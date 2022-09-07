function output = coinCellSectorGrid(params)

    % Components in order cathode on top, anode on bottom
    components = {'PositiveCurrentCollector', ...
                  'PositiveActiveMaterial', ...
                  'ElectrolyteSeparator', ...
                  'NegativeActiveMaterial', ...
                  'NegativeCurrentCollector'};

    for k = 1:numel(components)
        c = components{k};
        thickness(k) = params.thickness(c);
        diameter(k) = params.diameter(c);
        numCellLayers(k) = params.numCellLayers(c);
    end

    % Create tensor grid (extra large first radial element for slicing)
    nr = params.nR;
    R = max(diameter);
    r = linspace(0, R, nr);
    r(1) = -R;

    dz = thickness ./ numCellLayers;
    dz = rldecode(dz.', numCellLayers.');
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

    % Tag
    comptag = (1 : numel(components));
    tagdict = containers.Map(components, comptag);

    tag = -1*ones(G.cells.num, 1);
    zc = G.cells.centroids(:, 3);
    z = [-1 cumsum(thickness)];
    for k = 1:numel(z)-1
        t0 = z(k);
        t1 = z(k+1);
        idx = zc > t0 & zc < t1;
        c = components{k};
        tag(idx) = tagdict(c);
    end

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
    assert(numel(negativeExtCurrentFaces) == G.cartDims(1));
    assert(numel(positiveExtCurrentFaces) == G.cartDims(1));

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
