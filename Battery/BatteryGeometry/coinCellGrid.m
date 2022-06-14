function output = coinCellGrid(params)

    mrstModule add upr
    
    % Components in order cathode on top, anode on bottom
    components = {'PositiveCurrentCollector', ...
                  'PositiveActiveMaterial', ...
                  'ElectrolyteSeparator', ...
                  'NegativeActiveMaterial', ...
                  'NegativeCurrentCollector'};

    % Find diameter
    diameter = zeros(numel(components), 1);

    for k = 1:numel(components)
        key = components{k};
        diameter(k) = params.diameter(key);
    end

    R = 0.5 * max(diameter);
    bbox = [-1.2*R, -1.2*R; 1.2*R, -1.2*R; 1.2*R, 1.2*R; -1.2*R, 1.2*R];
    xc0 = mean(bbox);
    
    meshSize = params.meshSize;
    N = ceil(2 * pi * R / meshSize);

    outer_circle = circle(R, xc0, N);
    fc = outer_circle;

    % Find centers
    xc = repmat(xc0, numel(components), 1);
    if ~isempty(params.offset) & norm(params.offset) > 0
        for k = 1:numel(components)
            key = components{k};
            % Offset
            if (strcmpi(key, 'PositiveActiveMaterial') | strcmpi(key, 'NegativeActiveMaterial')) & ...
                    (~isempty(params.offset) & norm(params.offset) > 0)
                sign = 1;
                if strcmpi(key, 'NegativeActiveMaterial')
                    sign = -1;
                end
                xc(k, :) = xc0 + 0.5 * sign * params.offset;
            end
        end
    end

    % Create circles
    for k = 1:numel(components)
        key = components{k};
        if ismember(key, {'PositiveActiveMaterial', 'ElectrolyteSeparator', 'NegativeActiveMaterial'})
            d = params.diameter(key);
            r = 0.5 * d;
            n = ceil(2 * pi * r / meshSize);
            circ = circle(r, xc(k, :), n);
            fc = [fc, circ];
        end
    end

    G = pebiGrid2D(meshSize, max(bbox), 'faceConstraints', fc, 'polyBdr', bbox);
    G.faces = rmfield(G.faces, 'tag');

    %% Remove outside
    G = computeGeometry(G);
    d = vecnorm(G.cells.centroids - xc0, 2, 2);
    is_outside = d > R;
    G = removeCells(G, is_outside);

    %% Extrude
    for k = 1:numel(components)
        key = components{k};
        thickness(k) = params.thickness(key);
        numCellLayers(k) = params.numCellLayers(key);
    end
    dz = thickness ./ numCellLayers;
    dz = rldecode(dz.', numCellLayers.');
    G = makeLayeredGrid(G, dz);
    G = computeGeometry(G);

    % TODO: Chamfer negative cc
    %{

      |-------------|
      |             |
      |           _/
      |__________/

    %}

    % Tag
    comptag = 1 : numel(components);
    tagdict = containers.Map(components, comptag);

    thickness_accum = [0, cumsum(thickness)];
    zc = G.cells.centroids(:, 3);
    tag = -1*ones(G.cells.num, 1);
    
    for k = 1:numel(components)
        key = components{k};

        t0 = thickness_accum(k);
        t1 = thickness_accum(k+1);
        zidx = zc > t0 & zc < t1;
        
        rc = vecnorm(G.cells.centroids(:, 1:2) - xc(k), 2, 2);
        ridx = rc < 0.5 * diameter(k);
        
        idx = zidx & ridx;
        tag(idx) = tagdict(key);
    end
    
    %plotCellData(G, tag, G.cells.centroids(:,1)>0), view(3)    
    %keyboard;

    % Interfaces
    [bf, bc] = boundaryFaces(G);
    internal = (1:G.faces.num)';
    internal = setdiff(internal, bf);
    c1 = G.faces.neighbors(internal, 1);
    c2 = G.faces.neighbors(internal, 2);
    t1 = tag(c1);
    t2 = tag(c2);
    negidx = t1 == tagdict('NegativeActiveMaterial') & t2 == tagdict('NegativeCurrentCollector') | ...
             t2 == tagdict('NegativeActiveMaterial') & t1 == tagdict('NegativeCurrentCollector');
    negativeExtCurrentFaces = internal(find(negidx));
    posidx = t1 == tagdict('PositiveActiveMaterial') & t2 == tagdict('PositiveCurrentCollector') | ...
             t2 == tagdict('PositiveActiveMaterial') & t1 == tagdict('PositiveCurrentCollector');
    positiveExtCurrentFaces = internal(find(posidx));
    % assert(numel(negativeExtCurrentFaces) == G.cartDims(1));
    % assert(numel(positiveExtCurrentFaces) == G.cartDims(1));

    
    % Thermal cooling faces on all sides
    thermalCoolingFaces = bf;
    
    % top = G.faces.centroids(bf, 3) < 100*eps;
    % bottom = G.faces.centroids(bf, 3) > max(z) - 100*eps;
    % thermalCoolingFaces = [bf(top); bf(bottom)];
    %assert(numel(thermalCoolingFaces) == G.cartDims(1)*2);
    
    % No thermal exchange faces 
    % thermalExchangeFaces = find(G.faces.centroids(bf, 2) < 100*eps);
    % %assert(numel(thermalExchangeFaces) == G.cartDims(1)*G.cartDims(3));
    % thermalExchangeFacesTag(thermalExchangeFaces) = 1;
    
    % % Other side is found using looking at the correctly oriented
    % % normal direction
    % cc = G.cells.centroids(bc, :);
    % fc = G.faces.centroids(bf, :);
    % n = G.faces.normals(bf, :);
    % dotp = dot(fc - cc, n, 2);
    % neg = dotp < 0;
    % G.faces.normals(bf(neg),:) = -n(neg, :);
    % other = find(G.faces.normals(:,1) < 0);
    % %assert(numel(other) == numel(thermalExchangeFaces));
    % thermalExchangeFaces = [thermalExchangeFaces; other];
    % thermalExchangeFacesTag(other) = 2;

    thermalExchangeFaces = [];
    thermalExchangeFacesTag = [];


    output = params;

    output.G = G;
    output.tag = tag;
    output.tagdict = tagdict;

    output.negativeExtCurrentFaces = negativeExtCurrentFaces;
    output.positiveExtCurrentFaces = positiveExtCurrentFaces;

    % output.thermalCoolingFaces  = thermalCoolingFaces;
    % output.thermalExchangeFaces = thermalExchangeFaces;
    % output.thermalExchangeFacesTag = thermalExchangeFacesTag;

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
