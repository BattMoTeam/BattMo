function output = radialGrid(params)
    
    %% Parameters 
    % 
    % params structure with following fields
    % - nwindings : number of windings in the spiral
    % - r0        : "radius" at the middle
    % - widthDict : dictionary of widths for each component. The required key names for the dictionary are
    %                 - 'Separator'
    %                 - 'NegativeActiveMaterial'
    %                 - 'NegativeCurrentCollector'
    %                 - 'PositiveActiveMaterial'
    %                 - 'PositiveCurrentCollector'
    % - nrDict    : dicionary with number of cell in radial direction for each component (same keys as in widthDict).
    % - L         : length of the battery
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
    widthDict = params.widthDict ;
    nrDict    = params.nrDict;
    nas       = params.nas;
    L         = params.L;
    nL        = params.nL;


    %% component names
    
    compnames = {'PositiveActiveMaterial1', ...
                 'PositiveCurrentCollector', ...
                 'PositiveActiveMaterial2', ...
                 'ElectrolyteSeparator2', ...
                 'NegativeActiveMaterial2', ...
                 'NegativeCurrentCollector', ...
                 'NegativeActiveMaterial1', ...
                 'ElectrolyteSeparator1'};
    
    comptag = (1 : numel(compnames));
    tagdict = containers.Map(compnames, comptag);
    
    %% 
    widths = [widthDict('PositiveActiveMaterial'); ...
              widthDict('PositiveCurrentCollector'); ...
              widthDict('PositiveActiveMaterial'); ...
              widthDict('ElectrolyteSeparator'); ...
              widthDict('NegativeActiveMaterial'); ...
              widthDict('NegativeCurrentCollector'); ...
              widthDict('NegativeActiveMaterial'); ...
              widthDict('ElectrolyteSeparator')];
    
    nrs = [nrDict('PositiveActiveMaterial'); ...
           nrDict('PositiveCurrentCollector'); ...
           nrDict('PositiveActiveMaterial'); ...
           nrDict('ElectrolyteSeparator'); ...
           nrDict('NegativeActiveMaterial'); ...
           nrDict('NegativeCurrentCollector'); ...
           nrDict('NegativeActiveMaterial'); ...
           nrDict('ElectrolyteSeparator')];



    %% Grid setup

    r = rldecode(widths./nrs, nrs);
    r = repmat(r, nwindings, 1);
    r = r0 + [0; cumsum(r)];
    th = linspace(0, 2*pi, nas + 1);
    th = th(1 : end - 1);
    
    nR = size(r, 1) - 1;
    [R, TH]  = meshgrid(r, th);
    [px, py] = pol2cart(TH(:), R(:));
    [G, tR]  = buildRadialGrid([px, py], nas, nR);

    tbls = setupSimpleTables(G);
    
    celltbl = tbls.celltbl;
    
    celltbl = celltbl.addInd('indR', rldecode((1 : nR)', nas*ones(nR, 1)));
    celltbl = celltbl.addInd('indA', repmat((1 : nas)', nR, 1));
    
    ncomp = numel(widths);
    comptag = rldecode((1 : ncomp)', nrs);
    comptag = repmat(comptag, [nwindings, 1]);

    comptagtbl.tag = comptag;
    comptagtbl.indR = (1 : (sum(nrs)*nwindings))';
    comptagtbl = IndexArray(comptagtbl);

    celltagtbl = crossIndexArray(celltbl, comptagtbl, {'indR'});
    celltagtbl = sortIndexArray(celltagtbl, {'cells', 'tag'});

    tag = celltagtbl.get('tag');

    % Extrude battery in z-direction
    zwidths = (L/nL)*ones(nL, 1);
    G = makeLayeredGrid(G, zwidths);
    G = computeGeometry(G);
    
    tag = repmat(tag, [nL, 1]);

    % setup the standard tables
    tbls = setupSimpleTables(G);
    cellfacetbl = tbls.cellfacetbl;
    
    clear extfacetbl
    extfacetbl.faces = find(any(G.faces.neighbors == 0, 2));
    extfacetbl = IndexArray(extfacetbl);
    extcellfacetbl = crossIndexArray(extfacetbl, cellfacetbl, {'faces'});
    
    %%  recover the external faces that are inside the spiral
    % we get them using the Cartesian indexing

    [indA, indR, indZ] = ind2sub([nas, nR, nL], (1 : G.cells.num)');
    
    clear celltbl
    celltbl.cells = (1 : G.cells.num)';
    celltbl.indA = indA;
    celltbl.indR = indR;
    celltbl.indZ = indZ;
    celltbl = IndexArray(celltbl);

    % We add vertical (1) and horizontal (2) direction index for the faces (see makeLayeredGrid for the setup)
    
    nf = G.faces.num;
    clear facetbl
    facetbl.faces = (1 : nf)';
    dir = 2*ones(nf, 1);
    dir(1 : (nL + 1)*nas*nR) = 1;
    facetbl.dir = dir;
    facetbl = IndexArray(facetbl);
    
    scelltbl.indA = (1 : nas)';
    scelltbl.indR = 1*ones(nas, 1);
    scelltbl.dir = 2*ones(nas, 1);
    scelltbl = IndexArray(scelltbl);
    
    scelltbl = crossIndexArray(celltbl, scelltbl, {'indA', 'indR'});
    scelltbl = projIndexArray(scelltbl, 'cells');
    extscellfacetbl = crossIndexArray(scelltbl, extcellfacetbl, {'cells'});
    sfacetbl = projIndexArray(extscellfacetbl, {'faces'});
    
    sfacetbl = sfacetbl.addInd('dir', 2*ones(sfacetbl.num, 1));
    sfacetbl = crossIndexArray(sfacetbl, facetbl, {'faces', 'dir'});
    
    sfaces = sfacetbl.get('faces');
    
    clear sfacetbl
    sfacetbl.faces = sfaces;
    sfacetbl = IndexArray(sfacetbl);
    
    map = TensorMap();
    map.fromTbl = extfacetbl;
    map.toTbl = sfacetbl;
    map.mergefds = {'faces'};
    
    ind = map.getDispatchInd();
    
    thermalExchangeFaces = extfacetbl.get('faces');
    thermalExchangeFaces(ind) = [];
    
    
    %% recover faces on top and bottom for the current collector
    % we could do that using cartesian indices (easier)

    ccnames = {'PositiveCurrentCollector', 'NegativeCurrentCollector'};

    for ind = 1 : numel(ccnames)

        clear cccelltbl
        cccelltbl.cells = find(tag == tagdict(ccnames{ind}));
        cccelltbl = IndexArray(cccelltbl);

        extcccellfacetbl = crossIndexArray(extcellfacetbl, cccelltbl, {'cells'});
        extccfacetbl = projIndexArray(extcccellfacetbl, {'faces'});

        ccextfaces = extccfacetbl.get('faces');

        normals = G.faces.normals(ccextfaces, :);
        sgn = ones(numel(ccextfaces), 1);
        sgn(G.faces.neighbors(ccextfaces, 1) == 0) = -1;
        areas = G.faces.areas(ccextfaces, :);
        nnormals = sgn.*normals./areas;

        scalprod = bsxfun(@times, [0, 0, 1], nnormals);
        scalprod = sum(scalprod, 2);

        switch ccnames{ind}
          case 'PositiveCurrentCollector'
            ccfaces{ind} = ccextfaces(scalprod > 0.9);
          case 'NegativeCurrentCollector'
            ccfaces{ind} = ccextfaces(scalprod < -0.9);
          otherwise
            error('name not recognized');
        end
    end

    positiveExtCurrentFaces = ccfaces{1};
    negativeExtCurrentFaces = ccfaces{2};
   
    %% 
    detailedtag = tag;
    detailedtagdict = tagdict;

    tag = nan(G.cells.num);
    tagdict = containers.Map(...
        {'PositiveActiveMaterial'  , ...
         'PositiveCurrentCollector', ...
         'ElectrolyteSeparator'    , ...
         'NegativeActiveMaterial'  , ...
         'NegativeCurrentCollector'}, [1 : 5]');
    
    mappings = {{'PositiveActiveMaterial', {'PositiveActiveMaterial1', 'PositiveActiveMaterial2'}}, ...
                {'NegativeActiveMaterial', {'NegativeActiveMaterial1', 'NegativeActiveMaterial2'}}, ...
                {'ElectrolyteSeparator', {'ElectrolyteSeparator1', 'ElectrolyteSeparator2'}}, ...
                {'PositiveCurrentCollector', {'PositiveCurrentCollector'}}, ...
                {'NegativeCurrentCollector', {'NegativeCurrentCollector'}}};
    
    for ind1 = 1 : numel(mappings)
        mapping = mappings{ind1};
        tagvalue1 = tagdict(mapping{1});
        for ind2 = 1 : numel(mapping{2})
            tagvalue2 = detailedtagdict(mapping{2}{ind2});
            tag(detailedtag == tagvalue2) = tagvalue1;
        end
    end
    
    %% setup output structure
    output = params;
    output.G                       = G;
    output.tag                     = tag;
    output.tagdict                 = tagdict;
    output.positiveExtCurrentFaces = positiveExtCurrentFaces;
    output.negativeExtCurrentFaces = negativeExtCurrentFaces;
    output.thermalExchangeFaces    = thermalExchangeFaces;   
    
 
end