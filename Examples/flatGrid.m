function output = flatGrid(params)
    
    %% Parameters 
    % 
    % params structure with following fields
    % - nwindings : number of windings in the spiral
    % - r0        : "radius" at the middle
    % - widths    : vector of widths for each component indexed in following order (see tagdict below in code)
    %                 - positive current collector
    %                 - positive active material
    %                 - electrolyte separator 
    %                 - negative active material
    %                 - negative current collector
    % - L         : length of the battery
    % - nrs       : number of cell in radial direction for each component (same ordering as above).
    % - nas       : number of cells in the angular direction
    % - nL        : number of discretization cells in the longitudonal
    %
    % RETURNS
    %
    % - G       : grid
    % - tag     : cell-valued vector giving component number (indexing is given by tagdict)
    % - tagdict : dictionary giving the component number
    
    nwindings = params.nwindings;
    r0        = params.r0;
    widths    = params.widths;
    nrs       = params.nrs;
    nas       = params.nas;
    L         = params.L;
    nL        = params.nL;

    %% component names
    compnames = {'PositiveCurrentCollector', ...
                 'PositiveActiveMaterial'  , ...
                 'ElectrolyteSeparator'    , ...
                 'NegativeActiveMaterial'  , ...
                 'NegativeCurrentCollector'};
    
    comptag = (1 : numel(compnames));
    tagdict = containers.Map(compnames, comptag);

    
    %% Grid setup

    layerwidth = sum(widths);

    w = widths./nrs;
    w = rldecode(w, nrs);

    w = repmat(w, [nwindings, 1]);
    w = [0; cumsum(w)];

    h = linspace(0, 2*pi*r0, nas*nwindings + 1);

    nperlayer = sum(nrs);

    G = tensorGrid(h, w);

    n = numel(h) - 1;
    m = numel(w) - 1;

    tbls = setupSimpleTables(G);

    ncomp = numel(widths);
    comptag = rldecode((1 : ncomp)', nrs);
    comptag = repmat(comptag, [nwindings, 1]);

    comptagtbl.tag = comptag;
    comptagtbl.indj = (1 : (sum(nrs)*nwindings))';
    comptagtbl = IndexArray(comptagtbl);

    celltbl.cells = (1 : G.cells.num)';
    celltbl.indi = repmat((1 : nas*nwindings)', [sum(nrs)*nwindings, 1]);
    celltbl.indj = rldecode((1 : sum(nrs)*nwindings)', nas*nwindings*ones(sum(nrs)*nwindings, 1));
    celltbl = IndexArray(celltbl);

    celltagtbl = crossIndexArray(celltbl, comptagtbl, {'indj'});
    celltagtbl = sortIndexArray(celltagtbl, {'cells', 'tag'});

    tag = celltagtbl.get('tag');

    % Extrude battery in z-direction
    zwidths = (L/nL)*ones(nL, 1);
    G = makeLayeredGrid(G, zwidths);
    G = computeGeometry(G);
    
    G.faces = rmfield(G.faces, 'tag');
    
    tag = repmat(tag, [nL, 1]);

    % setup the standard tables
    tbls = setupSimpleTables(G);
    cellfacetbl = tbls.cellfacetbl;
    
    clear extfacetbl
    extfacetbl.faces = find(any(G.faces.neighbors == 0, 2));
    extfacetbl = IndexArray(extfacetbl);
    extcellfacetbl = crossIndexArray(extfacetbl, cellfacetbl, {'faces'});
    
    thermalExchangeFaces = extfacetbl.get('faces');
    
    [indi, indj, indk] = ind2sub([n, m, nL], (1 : G.cells.num)');
    
    clear celltbl
    celltbl.cells = (1 : G.cells.num)';
    celltbl.indi = indi;
    celltbl.indj = indj;
    celltbl.indk = indk;
    celltbl = IndexArray(celltbl);

    % We add horizontal (1) and vertical (2) direction index for the faces (see makeLayeredGrid for the setup)
    
    nf = G.faces.num;
    clear facetbl
    facetbl.faces = (1 : nf)';
    dir = 2*ones(nf, 1);
    dir(1 : (nL + 1)*n*m) = 1;
    facetbl.dir = dir;
    facetbl = IndexArray(facetbl);

    ccnames = {'PositiveCurrentCollector', 'NegativeCurrentCollector'};

    for ind = 1 : numel(ccnames)

        clear cccelltbl
        cccelltbl.cells = find(tag == tagdict(ccnames{ind}));
        switch ccnames{ind}
          case 'PositiveCurrentCollector'
            cccelltbl.indk = nL;
          case 'NegativeCurrentCollector'
            cccelltbl.indk = 1;
          otherwise
            error('name not recognized');
        end
        cccelltbl = IndexArray(cccelltbl);

        cccelltbl = crossIndexArray(cccelltbl, celltbl, {'cells', 'indk'});
        extcccellfacetbl = crossIndexArray(extcellfacetbl, cccelltbl, {'cells'});
        ccfacetbl = projIndexArray(extcccellfacetbl, {'faces'});
        
        ccfacetbl = ccfacetbl.addInd('dir', 1*ones(ccfacetbl.num, 1));
        
        ccfacetbl = crossIndexArray(ccfacetbl, facetbl, {'faces', 'dir'});

        ccfaces{ind} = ccfacetbl.get('faces');

    end

    positiveExtCurrentFaces = ccfaces{1};
    negativeExtCurrentFaces = ccfaces{2};
    
    %% setup output structure
    output = params;
    output.G                       = G;
    output.tag                     = tag;
    output.tagdict                 = tagdict;
    output.positiveExtCurrentFaces = positiveExtCurrentFaces;
    output.negativeExtCurrentFaces = negativeExtCurrentFaces;
    output.thermalExchangeFaces    = thermalExchangeFaces;   
    
end