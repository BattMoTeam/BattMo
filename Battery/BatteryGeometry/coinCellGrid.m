function output = coinCellGrid(params)

    mrstModule add upr

    % Components in order from z=0 (top) to z=zmax (bottom) with the
    % surrounding electrolyte last
    components = {'NegativeCurrentCollector', ...
                  'NegativeActiveMaterial', ...
                  'ElectrolyteSeparator', ...
                  'PositiveActiveMaterial', ...
                  'PositiveCurrentCollector', ...
                  'Electrolyte'};

    % Find parameters common for 3d grid and sector grid
    R = 0.5 * max(cell2mat(params.diameter.values));
    bbox = [-1.2*R, -1.2*R; 1.2*R, -1.2*R; 1.2*R, 1.2*R; -1.2*R, 1.2*R];
    xc0 = [0, 0];

    % Order is important here
    for k = 1:numel(components) - 1
        key = components{k};
        thickness(k) = params.thickness(key);
        numCellLayers(k) = params.numCellLayers(key);
    end
    dz = thickness ./ numCellLayers;
    dz = rldecode(dz.', numCellLayers.');

    comptag = 1 : numel(components);
    tagdict = containers.Map(components, comptag);

    if params.use_sector
        output = coinCellSectorGrid(params, components, R, bbox, xc0, thickness, dz, tagdict);
    else
        output = coinCell3dGrid(params, components, R, bbox, xc0, thickness, dz, tagdict);
    end

end


function output = coinCell3dGrid(params, components, R, bbox, xc0, thickness, dz, tagdict)

% Create constraints for the different diameters
    meshSize = params.meshSize;
    diameter = params.diameter;
    offset = params.offset;

    num_points = @(r) ceil(2 * pi * r / meshSize);
    
    if all(offset == 0)
        diams = unique(cell2mat(diameter.values));
        fc = cell(numel(diams), 1);
        for k = 1:numel(diams)
            r = 0.5 * diams(k);
            n = num_points(r);
            fc{k} = circle(r, xc0, n);
        end
    else
        % Circle data is [center, diameter] for all components
        circdata = zeros(numel(components)-1, 3);
        xc = containers.Map();
        for k = 1:numel(components)-1
            key = components{k};
            circdata(k, :) = [xc0, diameter(key)];
            if strcmp(key, 'PositiveActiveMaterial')
                circdata(k, 1:2) = circdata(k, 1:2) + 0.5 * offset;
            elseif strcmp(key, 'NegativeActiveMaterial')
                circdata(k, 1:2) = circdata(k, 1:2) - 0.5 * offset;
            end
        end

        % Create face constraints out of unique circle data
        cd = unique(circdata, 'rows');
        fc = cell(size(cd, 1), 1);
        for k = 1:size(cd, 1)
            xc = cd(k, 1:2);
            r = cd(k, 3);
            n = num_points(r);
            fc{k} = circle(r, xc, n);
        end
        keyboard;
    end
    
    fcfactor = 0.5;
    G = pebiGrid2D(meshSize, max(bbox), 'faceConstraints', [fc{:}], 'polyBdr', bbox, 'FCFactor', fcfactor);

    % Remove markers for faces constituting the outer circle
    G.faces = rmfield(G.faces, 'tag');

    %% Remove outside
    G = computeGeometry(G);
    d = vecnorm(G.cells.centroids - xc0, 2, 2);
    is_outside = d > R;
    G = removeCells(G, is_outside);

    %% Extrude
    G = makeLayeredGrid(G, dz);
    G = computeGeometry(G);

    %% Tag
    thickness_accum = [0, cumsum(thickness)];
    zc = G.cells.centroids(:, 3);
    tag = tagdict('Electrolyte') * ones(G.cells.num, 1);

    for k = 1:numel(components) - 1
        key = components{k};

        t0 = thickness_accum(k);
        t1 = thickness_accum(k+1);
        zidx = zc > t0 & zc < t1;

        rc = vecnorm(G.cells.centroids(:, 1:2) - xc0, 2, 2);
        ridx = rc < 0.5 * diameter(key);
        %ridx = true;

        idx = zidx & ridx;
        tag(idx) = tagdict(key);
    end
    assert(all(tag > -1));

    % Electrode faces where current is applied?
    [bf, bc] = boundaryFaces(G);
    minz = min(G.nodes.coords(:, 3)) + 100*eps;
    maxz = max(G.nodes.coords(:, 3)) - 100*eps;
    fmin = G.faces.centroids(:, 3) <= minz;
    fmax = G.faces.centroids(:, 3) >= maxz;
    assert(sum(fmin) == sum(fmax));

    cn = find(tag == tagdict('NegativeCurrentCollector'), 1);
    cp = find(tag == tagdict('PositiveCurrentCollector'), 1);

    if abs(cn - minz) < abs(cp - minz)
        negativeExtCurrentFaces = fmin;
        positiveExtCurrentFaces = fmax;
    else
        negativeExtCurrentFaces = fmax;
        positiveExtCurrentFaces = fmin;
    end

    % Thermal cooling faces on all sides
    thermalCoolingFaces = bf;



    output = params;

    output.G = G;
    output.tag = tag;
    output.tagdict = tagdict;

    output.negativeExtCurrentFaces = negativeExtCurrentFaces;
    output.positiveExtCurrentFaces = positiveExtCurrentFaces;

end


function output = coinCellSectorGrid(params, components, R, bbox, xc0, thickness, dz, tagdict)

% Create tensor grid (extra large first radial element for slicing)
    meshSize = 2*params.meshSize;
    nr = ceil(round(2*pi*R / meshSize));
    r0 = linspace(0, R, nr);
    diams = unique(cell2mat(params.diameter.values));
    r0 = unique([r0, 0.5*diams])';
    dr = diff(r0);
    r0(dr < 1e-10) = [];

    % smooth r
    r1 = laplacian_smoothing(r0(2:end-1), [r0(1); r0(end)]);
    r = r0;
    r(2:end-1) = r1;

    z = [0; cumsum(dz)];
    G = tensorGrid(r, [0, 1], z);

    % Shear top y nodes to mimic circles
    idx = abs(G.nodes.coords(:, 2) - 1) < eps & G.nodes.coords(:, 1) > 0;
    v = params.angle;
    shear = G.nodes.coords(idx, 1) * sin(v);
    G.nodes.coords(idx, 1) = G.nodes.coords(idx, 1) - shear;
    G = computeGeometry(G);

    % TODO: Chamfer negative cc
    %{

      |-------------|
      |             |
      |           _/
      |__________/

    %}
    % Slice and remove
    n = [-sin(v), cos(v), 0];
    pt = zeros(1, 3);
    [G, gix] = sliceGrid(G, pt, 'normal', n);
    if mrstPlatform('octave')
        legacy = true;
    else
        legacy = false;
    end
    m = markCutGrids(G, gix.new.faces, 'legacy', legacy);
    G = removeCells(G, m == 1);

    %% Tag
    thickness_accum = [-1, cumsum(thickness)];
    zc = G.cells.centroids(:, 3);
    tag = tagdict('Electrolyte') * ones(G.cells.num, 1);

    for k = 1:numel(components)-1
        key = components{k};

        t0 = thickness_accum(k);
        t1 = thickness_accum(k+1);
        zidx = zc > t0 & zc < t1;

        rc = vecnorm(G.cells.centroids(:, 1:2) - xc0, 2, 2);
        ridx = rc < 0.5 * params.diameter(key);

        idx = zidx & ridx;
        tag(idx) = tagdict(key);
    end

    % Electrode faces where current is applied?
    [bf, bc] = boundaryFaces(G);
    minz = min(G.nodes.coords(:, 3)) + 100*eps;
    maxz = max(G.nodes.coords(:, 3)) - 100*eps;
    fmin = G.faces.centroids(:, 3) <= minz;
    fmax = G.faces.centroids(:, 3) >= maxz;
    assert(sum(fmin) == sum(fmax));

    cn = find(tag == tagdict('NegativeCurrentCollector'), 1);
    cp = find(tag == tagdict('PositiveCurrentCollector'), 1);

    if abs(cn - minz) < abs(cp - minz)
        negativeExtCurrentFaces = fmin;
        positiveExtCurrentFaces = fmax;
    else
        negativeExtCurrentFaces = fmax;
        positiveExtCurrentFaces = fmin;
    end

    % Thermal cooling faces on top and bottom
    top = G.faces.centroids(bf, 3) < 100*eps;
    bottom = G.faces.centroids(bf, 3) > max(z) - 100*eps;
    thermalCoolingFaces = [bf(top); bf(bottom)];
    assert(numel(thermalCoolingFaces) == G.cartDims(1)*2);

    % Thermal exchange faces on the sides
    thermalExchangeFaces = find(G.faces.centroids(bf, 2) < 100*eps);
    assert(numel(thermalExchangeFaces) == G.cartDims(1)*G.cartDims(3));
    thermalExchangeFacesTag(thermalExchangeFaces) = 1;

    % Other side is found using looking at the correctly oriented
    % normal direction
    cc = G.cells.centroids(bc, :);
    fc = G.faces.centroids(bf, :);
    n = G.faces.normals(bf, :);
    dotp = dot(fc - cc, n, 2);
    neg = dotp < 0;
    G.faces.normals(bf(neg),:) = -n(neg, :);
    other = find(G.faces.normals(:,1) < 0);
    assert(numel(other) == numel(thermalExchangeFaces));
    thermalExchangeFaces = [thermalExchangeFaces; other];
    thermalExchangeFacesTag(other) = 2;



    output = params;

    output.G = G;
    output.tag = tag;
    output.tagdict = tagdict;

    output.negativeExtCurrentFaces = negativeExtCurrentFaces;
    output.positiveExtCurrentFaces = positiveExtCurrentFaces;

    output.thermalCoolingFaces  = thermalCoolingFaces;
    output.thermalExchangeFaces = thermalExchangeFaces;
    output.thermalExchangeFacesTag = thermalExchangeFacesTag;

end


function x = circle(r, xc, n)

    alpha = linspace(0, 2*pi, n + 1).';
    a = alpha(1:n);
    b = alpha(2:n+1);
    pts = @(v) [r*cos(v) + xc(1), r*sin(v) + xc(2)];
    xa = pts(a);
    xb = pts(b);
    for k = 1:n
        x{k} = [xa(k, :); xb(k, :)];
    end

end

function [y, A, rfix] = laplacian_smoothing(x, xfix, A, rfix)

% Return only smoothed y

    xx = [xfix; x];

    if nargin == 2
        nfix = (1:numel(xfix))';
        [xx, ii] = sort(xx);
        n = numel(xx);
        rfix = ismember(ii, nfix);

        A = diag(ones(n - 1, 1), 1);
        A = A + diag(ones(n - 1, 1),- 1);
        A(rfix, :) = 0;
        A = A + 2*diag(rfix);
        A = 0.5*A;
    end

    y = A * xx;
    y = y(~rfix);
end