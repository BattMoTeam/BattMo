function output = coinCellGrid(params)

    mrstModule add upr

    % Components in order from z=0 (top) to z=zmax (bottom)
    components = {'NegativeCurrentCollector', ...
                  'NegativeActiveMaterial', ...
                  'ElectrolyteSeparator', ...
                  'PositiveActiveMaterial', ...
                  'PositiveCurrentCollector'};

    R = 0.5 * max(cell2mat(params.diameter.values));
    bbox = [-1.2*R, -1.2*R; 1.2*R, -1.2*R; 1.2*R, 1.2*R; -1.2*R, 1.2*R];
    xc0 = mean(bbox);
    
    meshSize = params.meshSize;
    N = ceil(2 * pi * R / meshSize);

    outer_circle = circle(R, xc0, N);
    fc = outer_circle;

    fcfactor = 0.5;
    G = pebiGrid2D(meshSize, max(bbox), 'faceConstraints', fc, 'polyBdr', bbox, 'FCFactor', fcfactor);

    % Remove markers for faces constituting the outer circle
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
    
    %% Tag
    comptag = 1 : numel(components);
    tagdict = containers.Map(components, comptag);

    thickness_accum = [0, cumsum(thickness)];
    zc = G.cells.centroids(:, 3);
    %tag = tagdict('Electrolyte') * ones(G.cells.num, 1);
    tag = -1 * ones(G.cells.num, 1);
    
    for k = 1:numel(components)
        key = components{k};

        t0 = thickness_accum(k);
        t1 = thickness_accum(k+1);
        zidx = zc > t0 & zc < t1;

        %xc(k,:) = xc0;
        %rc = vecnorm(G.cells.centroids(:, 1:2) - xc(k), 2, 2);
        %ridx = rc < 0.5 * params.diameter(key);
        ridx = true;
        
        idx = zidx & ridx;
        tag(idx) = tagdict(key);
    end
    assert(all(tag > -1));
    
    % plotCellData(G, tag, G.cells.centroids(:,1)>0), view(3)    
    % keyboard;

    % % Interfaces
    % [bf, bc] = boundaryFaces(G);
    % internal = (1:G.faces.num)';
    % internal = setdiff(internal, bf);
    % c1 = G.faces.neighbors(internal, 1);
    % c2 = G.faces.neighbors(internal, 2);
    % t1 = tag(c1);
    % t2 = tag(c2);
    
    % negidx = t1 == tagdict('NegativeActiveMaterial') & t2 == tagdict('NegativeCurrentCollector') | ...
    %          t2 == tagdict('NegativeActiveMaterial') & t1 == tagdict('NegativeCurrentCollector');
    % negativeExtCurrentFaces = internal(find(negidx));

    % posidx = t1 == tagdict('PositiveActiveMaterial') & t2 == tagdict('PositiveCurrentCollector') | ...
    %          t2 == tagdict('PositiveActiveMaterial') & t1 == tagdict('PositiveCurrentCollector');
    % positiveExtCurrentFaces = internal(find(posidx));

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
    
    % cellneg = find(tag == tagdict('NegativeCurrentCollector'));
    % cellpos = find(tag == tagdict('PositiveCurrentCollector'));
    % c1 = G.faces.neighbors(bf, 1);
    % c1 = c1(c1>0);
    % c2 = G.faces.neighbors(bf, 2);
    % c2 = c2(c2>0);
    
    % % See if pos or neg is max z
    % zn = max(G.cells.centroids(cellneg, 3));
    % zp = max(G.cells.centroids(cellpos, 3));
    % if zn > zp
    %     maxz = G.faces.centroids(bf,3) >= zn;
    %     minz = G.faces.centroids(bf,3) <= min(G.cells.centroids(cellpos, 3));
    % end
    
    % negidx = tag(bc) == tagdict('NegativeCurrentCollector');
    
    % cn = G.faces.neighbors(bf, 1);
    % negidx = tag(cn) == tagdict('NegativeCurrentCollector');
    
    % keyboard;
    
    
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
