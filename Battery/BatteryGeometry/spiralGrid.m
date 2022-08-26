function output = spiralGrid(params)
    
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
    
    nwindings = params.nwindings;
    rInner        = params.rInner;
    widthDict = params.widthDict ;
    nrDict    = params.nrDict;
    nas       = params.nas;
    L         = params.L;
    nL        = params.nL;
    unifang   = params.angleuniform;
    tabparams = params.tabparams;

    
    tabcase = tabparams.tabcase;
    
    usetab = false;
    if ~strcmp('no tab', tabcase)
        usetab = true;
    end
        
    
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

    layerwidth = sum(widths);

    w = widths./nrs;
    w = rldecode(w, nrs);

    w = repmat(w, [nwindings, 1]);
    w = [0; cumsum(w)];

    if (unifang)
        h = linspace(0, 2*pi*rInner, nas + 1);
    else
        n1 = 15;
        n2 = nas - 2*n1;
        dtheta1 = 2*pi*0.02;
        dtheta2 = 2*pi - 2*dtheta1;
        dtheta = [repmat(dtheta1/n1, n1, 1); repmat(dtheta2/n2, n2, 1);repmat(dtheta1/n1, n1, 1)];
        theta = [0; cumsum(dtheta)];
        h = theta*rInner;
    end
    nperlayer = sum(nrs);

    cartG = tensorGrid(h, w);

    n = numel(h) - 1;
    m = numel(w) - 1;

    % shorcuts
    nx = nas;
    ny = sum(nrs);
    
    celltbl.cells = (1 : cartG.cells.num)';
    celltbl.indi = repmat((1 : nas)', [sum(nrs)*nwindings, 1]);
    celltbl.indj = rldecode((1 : sum(nrs)*nwindings)', nas*ones(sum(nrs)*nwindings, 1));
    celltbl.curvindi = celltbl.indi + nas*floor((celltbl.indj - 1)./ny);
    celltbl.curvindj = repmat(rldecode((1 : ny)', nx*ones(ny, 1)), nwindings, 1);
    celltbl = IndexArray(celltbl);
    
    
    % We roll the domain into a spirale
    x = cartG.nodes.coords(:, 1);
    y = cartG.nodes.coords(:, 2);

    theta = x./rInner;
    
    cartG.nodes.coords(:, 1) = (rInner + y + (theta/(2*pi))*layerwidth).*cos(theta);
    cartG.nodes.coords(:, 2) = (rInner + y + (theta/(2*pi))*layerwidth).*sin(theta);

    tbls = setupSimpleTables(cartG);

    % We add cartesian indexing for the nodes
    nodetbl.nodes = (1 : cartG.nodes.num)';
    nodetbl.indi = repmat((1 : (n + 1))', m + 1, 1);
    nodetbl.indj = rldecode((1 : (m + 1))', (n + 1)*ones(m + 1, 1));
    nodetbl = IndexArray(nodetbl);

    % We add cartesian indexing for the vertical faces (in original cartesian block)
    vertfacetbl.faces = (1 : (n + 1)*m)';
    vertfacetbl.indi = repmat((1 : (n + 1))', m, 1);
    vertfacetbl.indj = rldecode((1 : m)', (n + 1)*ones(m, 1));
    vertfacetbl = IndexArray(vertfacetbl);

    % Add structure to merge the nodes
    node2tbl.indi1 = ones(m - nperlayer + 1, 1);
    node2tbl.indj1 = ((nperlayer + 1) : (m + 1))';
    node2tbl.indi2 = (n + 1)*ones(m - nperlayer + 1, 1);
    node2tbl.indj2 = (1 : (m - nperlayer + 1))';
    node2tbl = IndexArray(node2tbl);

    gen = CrossIndexArrayGenerator();
    gen.tbl1 = nodetbl;
    gen.tbl2 = node2tbl;
    gen.replacefds1 = {{'indi', 'indi1'}, {'indj', 'indj1'}, {'nodes', 'nodes1'}};
    gen.mergefds = {'indi1', 'indj1'};

    node2tbl = gen.eval();

    gen = CrossIndexArrayGenerator();
    gen.tbl1 = nodetbl;
    gen.tbl2 = node2tbl;
    gen.replacefds1 = {{'indi', 'indi2'}, {'indj', 'indj2'}, {'nodes', 'nodes2'}};
    gen.mergefds = {'indi2', 'indj2'};

    node2tbl = gen.eval();

    node2tbl = sortIndexArray(node2tbl, {'nodes1', 'nodes2'});

    % Add structure to merge the faces
    face2tbl.indi1 = ones(m - nperlayer, 1);
    face2tbl.indj1 = nperlayer + (1 : (m - nperlayer))';
    face2tbl.indi2 = (n + 1)*ones(m - nperlayer, 1);
    face2tbl.indj2 = (1 : (m - nperlayer))';
    face2tbl = IndexArray(face2tbl);

    gen = CrossIndexArrayGenerator();
    gen.tbl1 = vertfacetbl;
    gen.tbl2 = face2tbl;
    gen.replacefds1 = {{'indi', 'indi1'}, {'indj', 'indj1'}, {'faces', 'faces1'}};
    gen.mergefds = {'indi1', 'indj1'};

    face2tbl = gen.eval();

    gen = CrossIndexArrayGenerator();
    gen.tbl1 = vertfacetbl;
    gen.tbl2 = face2tbl;
    gen.replacefds1 = {{'indi', 'indi2'}, {'indj', 'indj2'}, {'faces', 'faces2'}};
    gen.mergefds = {'indi2', 'indj2'};

    face2tbl = gen.eval();


    %% We setup the new indexing for the nodes

    nodetoremove = node2tbl.get('nodes2');
    newnodes = (1 : cartG.nodes.num)';
    newnodes(nodetoremove) = [];

    newnodetbl.newnodes = (1 : numel(newnodes))';
    newnodetbl.nodes = newnodes;
    newnodetbl = IndexArray(newnodetbl);

    gen = CrossIndexArrayGenerator();
    gen.tbl1 = node2tbl;
    gen.tbl2 = newnodetbl;
    gen.replacefds2 = {{'nodes', 'nodes1'}};
    gen.mergefds = {'nodes1'};

    node2tbl = gen.eval();

    newnodes = [newnodetbl.get('newnodes'); node2tbl.get('newnodes')];
    nodes = [newnodetbl.get('nodes'); node2tbl.get('nodes2')];

    clear newnodetbl;
    newnodetbl.newnodes = newnodes;
    newnodetbl.nodes = nodes;
    newnodetbl = IndexArray(newnodetbl);

    %% We setup the new indexing for the faces

    facetoremove = face2tbl.get('faces2');
    newfaces = (1 : cartG.faces.num)';
    newfaces(facetoremove) = [];

    clear facetbl
    newfacetbl.newfaces = (1 : numel(newfaces))';
    newfacetbl.faces = newfaces;
    newfacetbl = IndexArray(newfacetbl);

    gen = CrossIndexArrayGenerator();
    gen.tbl1 = face2tbl;
    gen.tbl2 = newfacetbl;
    gen.replacefds2 = {{'faces', 'faces1'}};
    gen.mergefds = {'faces1'};

    face2tbl = gen.eval();

    newfaces = [newfacetbl.get('newfaces'); face2tbl.get('newfaces')];
    faces = [newfacetbl.get('faces'); face2tbl.get('faces2')];

    allnewfacetbl.newfaces = newfaces;
    allnewfacetbl.faces = faces;
    allnewfacetbl = IndexArray(allnewfacetbl);

    %% we maps from old to new

    cellfacetbl = tbls.cellfacetbl;
    % we store the order previous to mapping. Here we just assumed that the grid is cartesian for simplicity
    cellfacetbl = cellfacetbl.addInd('order', repmat((1 : 4)', cartG.cells.num, 1));

    cellfacetbl = crossIndexArray(cellfacetbl, allnewfacetbl, {'faces'});
    cellfacetbl = sortIndexArray(cellfacetbl, {'cells', 'order' 'newfaces'});
    cellfacetbl = replacefield(cellfacetbl, {{'newfaces', 'faces'}});

    facenodetbl = tbls.facenodetbl;
    % facenodetbl = facenodetbl.addInd('order', repmat((1 : 2)', cartG.faces.num, 1));
    facenodetbl = crossIndexArray(facenodetbl, newfacetbl, {'faces'});
    facenodetbl = crossIndexArray(facenodetbl, newnodetbl, {'nodes'});

    facenodetbl = sortIndexArray(facenodetbl, {'newfaces',  'newnodes'});
    facenodetbl = replacefield(facenodetbl, {{'newfaces', 'faces'}, {'newnodes', 'nodes'}});

    clear nodes
    nodes.coords = cartG.nodes.coords(newnodetbl.get('nodes'), :);
    nodes.num = size(nodes.coords, 1);

    clear faces
    [~, ind] = rlencode(facenodetbl.get('faces'));
    faces.nodePos = [1; 1 + cumsum(ind)];
    faces.nodes = facenodetbl.get('nodes');
    faces.num = newfacetbl.num;
    faces.neighbors = []; % to avoid warning in computeGeometry
    
    clear cells
    [~, ind] = rlencode(cellfacetbl.get('cells'));
    cells.facePos = [1; 1 + cumsum(ind)];
    cells.faces = cellfacetbl.get('faces');
    cells.num = cartG.cells.num;

    G.cells = cells;
    G.faces = faces;
    G.nodes = nodes;
    G.griddim = 2;
    G.type = {'spiralGrid'};
    G = computeGeometry(G, 'findNeighbors', true);

    % plotGrid(G)

    ncomp = numel(widths);
    comptag = rldecode((1 : ncomp)', nrs);
    comptag = repmat(comptag, [nwindings, 1]);

    comptagtbl.tag = comptag;
    comptagtbl.indj = (1 : (sum(nrs)*nwindings))';
    comptagtbl = IndexArray(comptagtbl);

    celltagtbl = crossIndexArray(celltbl, comptagtbl, {'indj'});
    celltagtbl = sortIndexArray(celltagtbl, {'cells', 'tag'});

    tag = celltagtbl.get('tag');

    if usetab 
        
        w = widths./nrs;
        w = rldecode(w, nrs);
        indj0 = floor(nrDict('PositiveActiveMaterial') + nrDict('PositiveCurrentCollector')/2);
        cclinewidth = w(indj0);        
        
        cccelltbl.tag = tagdict('PositiveCurrentCollector');
        cccelltbl = IndexArray(cccelltbl);
        cccelltbl = crossIndexArray(cccelltbl, celltagtbl, {'tag'});
        cccelltbl = crossIndexArray(cccelltbl, celltbl, {'cells'});
        
        
        switch tabparams.tabcase

          case 'aligned tabs'
            
            windingnumbers = computeWindingNumbers(tabparams.fractions, rInner, layerwidth, nwindings);

            for ind = 1 :  numel(windingnumbers)
                
                indjwinding = indj0 + sum(nrs)*(windingnumbers(ind) - 1);
                clear cclinecelltbl
                cclinecelltbl.indj = repmat(indjwinding, nas, 1);
                cclinecelltbl.indi = (1 : nas)';
                cclinecelltbl = IndexArray(cclinecelltbl);
                
                cclinecelltbl = crossIndexArray(cclinecelltbl, celltbl, {'indi', 'indj'});
                cclinecelltbl = sortIndexArray(cclinecelltbl,  {'curvindi', 'cells'});            

                c = cclinecelltbl.get('cells');
                o = cclinecelltbl.get('curvindi');
                l = cumsum(G.cells.volumes(c)/cclinewidth);
                
                indl = find(l > tabparams.width, 1, 'first');
                tabwidths(ind) = l(indl);
                
                clear tabcelltbl
                tabcelltbl.curvindi = o(1 : indl);
                tabcelltbl = IndexArray(tabcelltbl);
                tabcelltbl = crossIndexArray(cccelltbl, tabcelltbl, {'curvindi'});
                tabcelltbl = projIndexArray(tabcelltbl, {'indi', 'indj'});
                
                tabcelltbls{ind} = tabcelltbl;
                
            end
            
            clear tabcelltbl
            tabcelltbl = tabcelltbls{1};
            for ind = 2 : numel(tabcelltbls)
                tabcelltbl = concatIndexArray(tabcelltbl, tabcelltbls{ind}, {'indi', 'indj'}, 'checkUnique', true);
            end
            
          case '1 tab'
            
            nmw = computeMiddleWinding(rInner, layerwidth, nwindings);
            indj0 = indj0 + sum(nrs)*(nmw - 1);
            cclinecelltbl.indj = repmat(indj0, nas, 1);
            cclinecelltbl.indi = (1 : nas)';
            cclinecelltbl = IndexArray(cclinecelltbl);
            
            cclinecelltbl = crossIndexArray(cclinecelltbl, celltbl, {'indi', 'indj'});
            cclinecelltbl = sortIndexArray(cclinecelltbl,  {'curvindi', 'cells'});            

            c = cclinecelltbl.get('cells');
            o = cclinecelltbl.get('curvindi');
            l = cumsum(G.cells.volumes(c)/cclinewidth);
            
            ind = find(l > tabparams.width, 1, 'first');
            tabwidths = l(ind);
            tabcelltbl.curvindi = o(1 : ind);
            tabcelltbl = IndexArray(tabcelltbl);
            
            tabcelltbl = crossIndexArray(cccelltbl, tabcelltbl, {'curvindi'});
            tabcelltbl = projIndexArray(tabcelltbl, {'indi', 'indj'});            

          case '3 tabs'

            indj0 = floor(nrDict('PositiveActiveMaterial') + nrDict('PositiveCurrentCollector')/2);
            
            cclinecelltbl.curvindj = indj0;
            cclinecelltbl = IndexArray(cclinecelltbl);
            cclinecelltbl = crossIndexArray(cclinecelltbl, celltbl, {'curvindj'});
            cclinecelltbl = crossIndexArray(cclinecelltbl, celltbl, {'indi', 'indj', 'curvindi'});
            cclinecelltbl = sortIndexArray(cclinecelltbl,  {'curvindi', 'cells'});


            
            c = cclinecelltbl.get('cells');
            o = cclinecelltbl.get('curvindi');
            l = cumsum(G.cells.volumes(c)/cclinewidth);

            tabwidths = tabparams.widths;
            ntab = numel(tabwidths);
            assert(ntab == 3); % only this is supported. It could be easily changed
            
            tabcells = [];
            
            for tabind = 1 : ntab
                
                switch tabind
                  case 1
                    tabcenter = pi*rInner;
                    ind = (l > tabcenter - tabwidths(tabind)/2);
                    ind = ind & (l < tabcenter + tabwidths(tabind)/2);
                  case 2
                    tabcenter = l(end)/2;
                    ind = (l > tabcenter - tabwidths(tabind)/2);
                    ind = ind & (l < tabcenter + tabwidths(tabind)/2);
                  case 3
                    tabcenter = l(end);
                    ind = (l > l(end) - tabwidths(tabind));
                end
                
                if isempty(ind)
                    ind = find(l >= tabcenter, 1, 'first');
                    if ind > 0 && abs(l(ind - 1) - tabcenter)  < abs(l(ind) - tabcenter);
                        ind = ind - 1
                    end
                end
                
                otab = o(ind);
                
                clear selcurvindtbl;
                selcurvindtbl.curvindi = otab;
                selcurvindtbl = IndexArray(selcurvindtbl);
                
                tabcelltbl = crossIndexArray(cccelltbl, selcurvindtbl, {'curvindi'});
                c = tabcelltbl.get('cells');
                
                tabcells = vertcat(tabcells, c);
            end
            
            clear tabcelltbl
            tabcelltbl.cells = tabcells;
            tabcelltbl = IndexArray(tabcelltbl);
            tabcelltbl = crossIndexArray(tabcelltbl, celltbl, {'cells'});
            tabcelltbl = projIndexArray(tabcelltbl, {'indi', 'indj'});

          otherwise
            error('tabcase not recognized');
        end    
    end
    
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
    % We get them using the Cartesian indexing

    cellindktbl.indk = (1 : nL)';
    cellindktbl = IndexArray(cellindktbl);
    celltbl = celltbl.removeInd({'cells'});
    celltbl = crossIndexArray(cellindktbl, celltbl, {});
    celltbl = sortIndexArray(celltbl, {'indk', 'indj', 'indi', 'curvindi', 'curvindj'});
    celltbl = celltbl.addInd('cells', (1 : G.cells.num)');

    % We add vertical (2) and horizontal (1) direction index for the faces (see makeLayeredGrid for the setup)
    
    nf = G.faces.num;
    clear facetbl
    facetbl.faces = (1 : nf)';
    dir = 2*ones(nf, 1);
    dir(1 : (nL + 1)*n*m) = 1;
    facetbl.dir = dir;
    facetbl = IndexArray(facetbl);
    
    %% We extract the faces at the exterior for thermal exchange, using Cartesian indexing
    
    scelltbl.indi = (1: n)';
    scelltbl.indj = 1*ones(n, 1);
    scelltbl.dir  = 2*ones(n, 1);
    scelltbl = IndexArray(scelltbl);
    
    scelltbl = crossIndexArray(celltbl, scelltbl, {'indi', 'indj'});
    scelltbl = projIndexArray(scelltbl, 'cells');
    extscellfacetbl = crossIndexArray(scelltbl, extcellfacetbl, {'cells'});
    sfacetbl = projIndexArray(extscellfacetbl, {'faces'});
    
    sfacetbl = sfacetbl.addInd('dir', 2*ones(sfacetbl.num, 1));
    sfacetbl = crossIndexArray(sfacetbl, facetbl, {'faces', 'dir'});
    
    sfaces = sfacetbl.get('faces'); % external vertical faces of the internal layer (inside the roll). It does not
                                    % include the border of the layer.
    
    clear scelltbl
    nnrs = sum(nrs);
    scelltbl.indi = ones(nnrs, 1);
    scelltbl.indj = (1 : nnrs)';
    scelltbl = IndexArray(scelltbl);
    
    scelltbl = crossIndexArray(celltbl, scelltbl, {'indi', 'indj'});
    scelltbl = projIndexArray(scelltbl, 'cells');
    extscellfacetbl = crossIndexArray(scelltbl, extcellfacetbl, {'cells'});
    sfacetbl = projIndexArray(extscellfacetbl, {'faces'});
    
    sfacetbl = sfacetbl.addInd('dir', 2*ones(sfacetbl.num, 1));
    sfacetbl = crossIndexArray(sfacetbl, facetbl, {'faces', 'dir'});

    sfaces = [sfaces; sfacetbl.get('faces')]; % We add the external faces of the border of the internal layer.
    
    % some faces have been taken twice
    sfaces = unique(sfaces);

    clear nonthfacetbl
    nonthfacetbl.faces = sfaces;
    nonthfacetbl = IndexArray(nonthfacetbl);
    
    map = TensorMap();
    map.fromTbl = extfacetbl;
    map.toTbl = nonthfacetbl;
    map.mergefds = {'faces'};
    
    ind = map.getDispatchInd();
    
    thermalExchangeFaces = extfacetbl.get('faces');
    thermalExchangeFaces(ind) = [];
    
    %% We capture the top faces of the internal 
    clear scelltbl
    topcelltbl.indk = nL;
    topcelltbl = IndexArray(topcelltbl);
    
    topcelltbl = crossIndexArray(topcelltbl, celltbl, {'indk'});
    topcellfacetbl = crossIndexArray(extcellfacetbl, topcelltbl, {'cells'});
    topfacetbl = projIndexArray(topcellfacetbl, {'faces'});
    
    hfacetbl.dir = 1;
    hfacetbl = IndexArray(hfacetbl);
    hfacetbl = crossIndexArray(facetbl, hfacetbl, {'dir'});
    
    topfacetbl = crossIndexArray(topfacetbl, hfacetbl, {'faces'});
    
    thfacetbl.faces = thermalExchangeFaces;
    thfacetbl = IndexArray(thfacetbl);
    
    map = TensorMap();
    map.fromTbl = thfacetbl;
    map.toTbl = topfacetbl;
    map.mergefds =  {'faces'};
    
    ind = map.getDispatchInd();
    
    % Thermal faces tag is equal to 1 if top faces and equal to 2 otherwise.
    thermalExchangeFacesTag = 2*ones(thfacetbl.num, 1);
    thermalExchangeFacesTag(ind) = 1;

    %% recover faces on top and bottom for the current collector

    ccnames = {'PositiveCurrentCollector', 'NegativeCurrentCollector'};

    for ind = 1 : numel(ccnames)

        ccname = ccnames{ind};
        clear cccelltbl
        cccelltbl.cells = find(tag == tagdict(ccnames{ind}));
        cccelltbl = IndexArray(cccelltbl);
        if strcmp(ccname, 'PositiveCurrentCollector')
            cccelltbl = cccelltbl.addInd('indk', nL*ones(cccelltbl.num, 1));
        else
            cccelltbl = cccelltbl.addInd('indk', 1*ones(cccelltbl.num, 1));
        end
        cccelltbl = crossIndexArray(cccelltbl, celltbl, {'cells', 'indk'});
        
        if usetab && strcmp(ccname, 'PositiveCurrentCollector') 
            cccelltbl = crossIndexArray(cccelltbl, tabcelltbl, {'indi', 'indj'});
        end
        ccfacetbl = crossIndexArray(cccelltbl, cellfacetbl, {'cells'});
        ccfacetbl = projIndexArray(ccfacetbl, {'faces'});
        ccfacetbl = ccfacetbl.addInd('dir', ones(ccfacetbl.num, 1));
        
        ccfacetbl = crossIndexArray(facetbl, ccfacetbl, {'faces', 'dir'});
        ccfacetbl = crossIndexArray(ccfacetbl, extfacetbl, {'faces'});
        
        ccfaces{ind} = ccfacetbl.get('faces');

    end

    positiveExtCurrentFaces = ccfaces{1};
    negativeExtCurrentFaces = ccfaces{2};
    
    %% 
    detailedtag = tag;
    detailedtagdict = tagdict;

    tag = nan(G.cells.num, 1);
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
    output.thermalExchangeFacesTag = thermalExchangeFacesTag;
    
    switch tabcase
      case '1 tab'
        output.tabwidths = tabwidths;
      case 'aligned tabs'
        output.tabwidths = tabwidths;
        output.windingnumbers = windingnumbers;
    end
    
    w = widths./nrs;
    w = rldecode(w, nrs);
    h = L/nL;
    h = repmat(h, nL, 1);
    
    output.widthLayer   = w;
    output.nWidthLayer  = nrs;
    output.heightLayer  = h;
    output.nHeightLayer = nL;
    output.celltbl      = celltbl;
    
end
    




%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
