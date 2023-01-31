function output = coinCellGrid(params)

    % Components in order from z=0 (top) to z=zmax (bottom) with the
    % surrounding electrolyte last
    compNames = {'NegativeCurrentCollector', ...
                 'NegativeActiveMaterial', ...
                 'ElectrolyteSeparator', ...
                 'PositiveActiveMaterial', ...
                 'PositiveCurrentCollector', ...
                 'Electrolyte'};
    assert(numel(intersect(compNames, params.compDims.Row)) == 5);
    
    thickness = params.compDims.thickness;
    numCellLayers = params.compDims.numCellLayers;
    dz = thickness ./ numCellLayers;

    % Repeat dz for each layer
    dz = rldecode(dz, numCellLayers);
    
    comptag = 1 : numel(compNames);
    tagdict = containers.Map(compNames, comptag);

    output = coinCell3dGrid(params, compNames, dz, tagdict);

end


function output = coinCell3dGrid(params, compNames, dz, tagdict)

    compDims = params.compDims;

    % Extract basic dimensions
    R = 0.5 * max(compDims.diameter);
    bbox = [-1.2*R, -1.2*R; 1.2*R, -1.2*R; 1.2*R, 1.2*R; -1.2*R, 1.2*R];
    xc0 = [0, 0];
    
    % Radial points
    hr = R / params.numRadial;
    diams = unique(compDims.diameter);
    r = 0.5 * diams;
    dr = diff(r);
    x1 = linspace(0, r(1), ceil(r(1)/hr));

    % Make center cell smaller:
    x1 = [0 0.25*x1(2) x1(2:end)];

    % Make sure we have 3 nodes.
    %nr = max(3, ceil(dr/hr));
    nr = ceil(dr/hr);

    for k = 2:numel(r)
        x1 = [x1, linspace(r(k-1), r(k), nr(k-1))];
    end
    x1 = unique(x1)';
    x1(:,2) = zeros(size(x1,1), 1);

    % Curves for transfinite interpolation
    x1(1,:) = [];
    s = (x1(:,1)-x1(1,1))/(x1(end,1)-x1(1,1));
    f1 = x1;
    f2 = f1;
    t = linspace(0, 1, params.numAngular+1)';
    c = [cos(t*2*pi), sin(t*2*pi)];
    g1 = f1(1,1)*c;
    g2 = f1(end,1)*c;

    G = transfiniteGrid(s, t, f1, f2, g1, g2, 'fill', true);

    % Extrude
    G = makeLayeredGrid(G, dz);
    G = computeGeometry(G);

    %% Tag
    thickness = compDims.thickness;
    thickness_accum = [0; cumsum(thickness)];
    zc = G.cells.centroids(:, 3);
    tag = tagdict('Electrolyte') * ones(G.cells.num, 1);

    for k = 1:numel(compNames) 
        key = compNames{k};
        if ~strcmp(key, 'Electrolyte')
            t0 = thickness_accum(k);
            t1 = thickness_accum(k+1);
            zidx = zc > t0 & zc < t1;
            
            rc = vecnorm(G.cells.centroids(:, 1:2) - xc0, 2, 2);
            ridx = rc < 0.5 * compDims{key, 'diameter'};
            
            idx = zidx & ridx;
            tag(idx) = tagdict(key);
        end
    end
    assert(all(tag > -1));
    
    %% Electrode faces where current is applied
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

    % Output
    output = params;
    output.G = G;
    output.tag = tag;
    output.tagdict = tagdict;
    output.negativeExtCurrentFaces = negativeExtCurrentFaces;
    output.positiveExtCurrentFaces = positiveExtCurrentFaces;

end
