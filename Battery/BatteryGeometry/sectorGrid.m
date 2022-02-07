function output = sectorGrid(params)
    
    %% Parameters 
    % 
    % params structure with following fields
    % - nwindings : number of windings in the spiral
    % - rInner        : "radius" at the middle
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
    

    output = radialGrid(params);
    
    
    nwindings = params.nwindings;
    nrDict    = params.nrDict;
    nas       = params.nas;
    nL        = params.nL;

    nrs = [nrDict('PositiveActiveMaterial'); ...
           nrDict('PositiveCurrentCollector'); ...
           nrDict('PositiveActiveMaterial'); ...
           nrDict('ElectrolyteSeparator'); ...
           nrDict('NegativeActiveMaterial'); ...
           nrDict('NegativeCurrentCollector'); ...
           nrDict('NegativeActiveMaterial'); ...
           nrDict('ElectrolyteSeparator')];
        
    nR = sum(nrs)*nwindings;
    
    radG                    = output.G;
    tag                     = output.tag;
    tagdict                 = output.tagdict;
    positiveExtCurrentFaces = output.positiveExtCurrentFaces;
    negativeExtCurrentFaces = output.negativeExtCurrentFaces;
    thermalExchangeFaces    = output.thermalExchangeFaces;
    
    [indA, indR, indZ] = ind2sub([nas, nR, nL], (1 : radG.cells.num)');
    celltbl.cells = (1 : radG.cells.num)';
    celltbl.indA = indA;
    celltbl.indR = indR;
    celltbl.indZ = indZ;
    celltbl = IndexArray(celltbl);
    
    scelltbl.indA = 1;
    scelltbl = IndexArray(scelltbl);
    
    scelltbl = crossIndexArray(scelltbl, celltbl, {'indA'});
    
    scells = scelltbl.get('cells');
    
    rcells = (1 : radG.cells.num)';
    rcells(scells) = [];
    
    [G, cellmap, facemap, nodemap] = removeCells(radG, rcells);
    
    tag = tag(cellmap);
    
    tbls = setupSimpleTables(radG);
    radfacetbl = tbls.facetbl;
    tbls = setupSimpleTables(G);
    facetbl = tbls.facetbl;
    facetbl = facetbl.addInd('radfaces', facemap);

    map = TensorMap();
    map.fromTbl = radfacetbl;
    map.toTbl = facetbl;
    map.replaceFromTblfds = {{'faces', 'radfaces'}};
    map.mergefds = {'radfaces'};
    map = map.setup();

    sfaces = zeros(radfacetbl.num, 1);
    sfaces(positiveExtCurrentFaces) = 1;
    sfaces = map.eval(sfaces);
    positiveExtCurrentFaces = find(sfaces);
    
    sfaces = zeros(radfacetbl.num, 1);
    sfaces(negativeExtCurrentFaces) = 1;
    sfaces = map.eval(sfaces);
    negativeExtCurrentFaces = find(sfaces);
    
    sfaces = zeros(radfacetbl.num, 1);
    sfaces(thermalExchangeFaces) = 1;
    sfaces = map.eval(sfaces);
    thermalExchangeFaces = find(sfaces);
    
    %% setup output structure
    output = params;
    output.G                       = G;
    output.tag                     = tag;
    output.tagdict                 = tagdict;
    output.positiveExtCurrentFaces = positiveExtCurrentFaces;
    output.negativeExtCurrentFaces = negativeExtCurrentFaces;
    output.thermalExchangeFaces    = thermalExchangeFaces;   
 
end