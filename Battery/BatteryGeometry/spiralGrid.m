function output = spiralGrid(params)

    %% Parameters
    %
    % params structure with following fields
    % - nwindings : number of windings in the spiral
    % - rInner        : "radius" at the middle
    % - widthDict : dictionary of widths for each component. The required key names for the dictionary are
    %                 - 'Separator'
    %                 - 'NegativeCoating'
    %                 - 'NegativeCurrentCollector'
    %                 - 'PositiveCoating'
    %                 - 'PositiveCurrentCollector'
    % - nrDict    : dicionary with number of cell in radial direction for each component (same keys as in widthDict).
    % - L         : length of the battery
    % - nas       : number of cells in the angular direction
    % - nL        : number of discretization cells in the longitudonal
    % - refLcoef  : coefficient use in refinement at top/bottom (if empty, no refinement)
    % - tabparams : structure with two fields
    %     - NegativeElectrode : for the negative current collector
    %     - PositiveElectrode : for the positive current collector
    %               Each of this structure contains the fields (see json schema battmoDir()/Utilities/JsonSchemas/Geometry.schema.json)
    %        - usetab    : boolean
    %        - fractions : length fraction where the tabs are located (one value per tab), only used if usetab is true
    %        - width     : width of the tabs (for the moment, only a single value), only used if usetab is true
    %
    % RETURNS
    %
    % - output structure with fields
    %    - G
    %    - tag
    %    - tagdict
    %    - positiveExtCurrentFaces
    %    - negativeExtCurrentFaces
    %    - thermalExchangeFaces
    %    - thermalExchangeFacesTag
    %    - widthLayer
    %    - nWidthLayer
    %    - heightLayer
    %    - nHeightLayer
    %    - celltbl

    nwindings  = params.nwindings;
    rInner     = params.rInner;
    widthDict  = params.widthDict ;
    nrDict     = params.nrDict;
    nas        = params.nas;
    L          = params.L;
    nL         = params.nL;
    unifang    = params.angleuniform;
    tabparams  = params.tabparams;

    exteriorNegativeElectrodeLayer = params.exteriorNegativeElectrodeLayer; % boolean

    
    %% component names

    compnames = {'PositiveCoating1'        , ...
                 'PositiveCurrentCollector', ...
                 'PositiveCoating2'        , ...
                 'Separator2'              , ...
                 'NegativeCoating2'        , ...
                 'NegativeCurrentCollector', ...
                 'NegativeCoating1'        , ...
                 'Separator1'};

    comptag = (1 : numel(compnames));
    tagdict = containers.Map(compnames, comptag);

    %%
    widths = [widthDict('PositiveCoating')          ; ...
              widthDict('PositiveCurrentCollector') ; ...
              widthDict('PositiveCoating')          ; ...
              widthDict('Separator')                ; ...
              widthDict('NegativeCoating')          ; ...
              widthDict('NegativeCurrentCollector') ; ...
              widthDict('NegativeCoating')          ; ...
              widthDict('Separator')];

    nrs = [nrDict('PositiveCoating')          ; ...
           nrDict('PositiveCurrentCollector') ; ...
           nrDict('PositiveCoating')          ; ...
           nrDict('Separator')                ; ...
           nrDict('NegativeCoating')          ; ...
           nrDict('NegativeCurrentCollector') ; ...
           nrDict('NegativeCoating')          ; ...
           nrDict('Separator')];

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

    celltbl.cells    = (1 : cartG.cells.num)';
    celltbl.indi     = repmat((1 : nas)', [sum(nrs)*nwindings, 1]);
    celltbl.indj     = rldecode((1 : sum(nrs)*nwindings)', nas*ones(sum(nrs)*nwindings, 1));
    celltbl.curvindi = celltbl.indi + nas*floor((celltbl.indj - 1)./ny);
    celltbl.curvindj = repmat(rldecode((1 : ny)', nx*ones(ny, 1)), nwindings, 1);
    celltbl = IndexArray(celltbl);

    % We roll the domain into a spirale
    x = cartG.nodes.coords(:, 1);
    y = cartG.nodes.coords(:, 2);

    theta = x./rInner;

    cartG.nodes.coords(:, 1) = (rInner + y + (theta/(2*pi))*layerwidth).*cos(theta);
    cartG.nodes.coords(:, 2) = (rInner + y + (theta/(2*pi))*layerwidth).*sin(theta);

    tbls = setupTables(cartG);

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

    if isfield(tabparams, 'PositiveElectrode') && tabparams.PositiveElectrode.usetab
        petab = tabparams.PositiveElectrode;
    else
        petab = [];
    end

    if isfield(tabparams, 'NegativeElectrode') && tabparams.NegativeElectrode.usetab
        netab = tabparams.NegativeElectrode;
    else
        netab = [];
    end

    if ~isempty(petab) || ~isempty(netab)
        clear data
        data.nas     = nas;
        data.nrs     = nrs;
        data.celltbl = celltbl;
        data.G       = G;
        data.tagdict = tagdict;
    end

    if ~isempty(petab)

        w = widths./nrs;
        w = rldecode(w, nrs);
        indj0 = floor(nrDict('PositiveCoating') + nrDict('PositiveCurrentCollector')/2);
        cclinewidth = w(indj0);

        clear cccelltbl;
        cccelltbl.tag = tagdict('PositiveCurrentCollector');
        cccelltbl = IndexArray(cccelltbl);
        cccelltbl = crossIndexArray(cccelltbl, celltagtbl, {'tag'});
        cccelltbl = crossIndexArray(cccelltbl, celltbl, {'cells'});

        fractions = petab.fractions;

        windingnumbers = computeWindingNumbers(fractions, rInner, layerwidth, nwindings);

        pe_tabcelltbl = getTabCellTbl(cccelltbl     , ...
                                      cclinewidth   , ...
                                      indj0         , ...
                                      petab         , ...
                                      windingnumbers, ...
                                      data);

    end

    if ~isempty(netab)

        w = widths./nrs;
        w = rldecode(w, nrs);
        indj0 = floor(nrDict('NegativeCoating') + nrDict('NegativeCurrentCollector')/2);
        cclinewidth = w(indj0);

        clear cccelltbl;
        cccelltbl.tag = tagdict('NegativeCurrentCollector');
        cccelltbl = IndexArray(cccelltbl);
        cccelltbl = crossIndexArray(cccelltbl, celltagtbl, {'tag'});
        cccelltbl = crossIndexArray(cccelltbl, celltbl, {'cells'});

        fractions = netab.fractions;

        windingnumbers = computeWindingNumbers(fractions, rInner, layerwidth, nwindings);

        ne_tabcelltbl = getTabCellTbl(cccelltbl     , ...
                                      cclinewidth   , ...
                                      indj0         , ...
                                      netab         , ...
                                      windingnumbers, ...
                                      data);

    end

    % Extrude battery in z-direction
    zwidths = (L/nL)*ones(nL, 1);

    % apply refinement if given
    if (isa(params, 'BatteryGenerator') | isfield(params, 'refLcoef')) && ~isempty(params.refLcoef)
        alpha = params.refLcoef;
        z = [0; cumsum(zwidths)];
        z = L*0.5*(atan(alpha*(z/L - 1/2))/atan(alpha/2) + 1);
        zwidths = diff(z);
    end
    
    G = makeLayeredGrid(G, zwidths);
    G = computeGeometry(G);

    tag = repmat(tag, [nL, 1]);

    % setup the standard tables
    tbls = setupTables(G);
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
    ccfaces = cell(numel(ccnames), 1);

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

        if ~isempty(petab) && strcmp(ccname, 'PositiveCurrentCollector')
            cccelltbl = crossIndexArray(cccelltbl, pe_tabcelltbl, {'indi', 'indj'});
        end

        if ~isempty(netab) && strcmp(ccname, 'NegativeCurrentCollector')
            cccelltbl = crossIndexArray(cccelltbl, ne_tabcelltbl, {'indi', 'indj'});
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
        {'PositiveCoating'         , ...
         'PositiveCurrentCollector', ...
         'Separator'               , ...
         'NegativeCoating'         , ...
         'NegativeCurrentCollector'}, (1 : 5)');

    mappings = {{'PositiveCoating'          , {'PositiveCoating1', 'PositiveCoating2'}}, ...
                {'NegativeCoating'          , {'NegativeCoating1', 'NegativeCoating2'}}, ...
                {'Separator'                , {'Separator1', 'Separator2'}}            , ...
                {'PositiveCurrentCollector' , {'PositiveCurrentCollector'}}            , ...
                {'NegativeCurrentCollector' , {'NegativeCurrentCollector'}}};

    for ind1 = 1 : numel(mappings)
        mapping = mappings{ind1};
        tagvalue1 = tagdict(mapping{1});
        for ind2 = 1 : numel(mapping{2})
            tagvalue2 = detailedtagdict(mapping{2}{ind2});
            tag(detailedtag == tagvalue2) = tagvalue1;
        end
    end

    if exteriorNegativeElectrodeLayer

        sel_indj_tbl.indj = ((sum(nrs(1 : 4)) + 1) : sum(nrs)*nwindings)';
        sel_indj_tbl = IndexArray(sel_indj_tbl);

        celltbl = crossIndexArray(celltbl, sel_indj_tbl, {'indj'});

        c = celltbl.get('cells');

        [G, cellmap, facemap, nodemap] =  extractSubgrid(G, c);

        % map cell
        oldnewcelltbl.oldcells = cellmap;
        oldnewcelltbl.newcells = (1 : G.cells.num)';
        oldnewcelltbl = IndexArray(oldnewcelltbl);

        gen = CrossIndexArrayGenerator();
        gen.tbl1        = celltbl;
        gen.tbl2        = oldnewcelltbl;
        gen.replacefds1 = {{'cells', 'oldcells'}};
        gen.replacefds2 = {{'newcells', 'cells'}};
        gen.mergefds    = {'oldcells'};

        celltbl = gen.eval();
        celltbl = celltbl.removeInd({'oldcells'});
        
        tag = tag(cellmap);

        % map faces
        newoldfacetbl.newfaces = (1 : G.faces.num)';
        newoldfacetbl.oldfaces = facemap;
        newoldfacetbl = IndexArray(newoldfacetbl);

        positiveExtCurrentFaces = mapToNewFaces(positiveExtCurrentFaces, newoldfacetbl);
        negativeExtCurrentFaces = mapToNewFaces(negativeExtCurrentFaces, newoldfacetbl);
        
        thexch_face_tbl.oldfaces = thermalExchangeFaces;
        thexch_face_tbl = IndexArray(thexch_face_tbl);

        [thexch_face_tbl, indstruct] = crossIndexArray(thexch_face_tbl, newoldfacetbl, {'oldfaces'});
        
        thermalExchangeFaces = thexch_face_tbl.get('newfaces');
        thermalExchangeFacesTag = thermalExchangeFacesTag(indstruct{1}.inds);
        
    end
    
    %% setup output structure
    output.G                       = G;
    output.tag                     = tag;
    output.tagdict                 = tagdict;
    output.positiveExtCurrentFaces = positiveExtCurrentFaces;
    output.negativeExtCurrentFaces = negativeExtCurrentFaces;
    output.thermalExchangeFaces    = thermalExchangeFaces;
    output.thermalExchangeFacesTag = thermalExchangeFacesTag;

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

function newfaces = mapToNewFaces(oldfaces, newoldfacetbl)
    
    selfacetbl.oldfaces = oldfaces;
    selfacetbl = IndexArray(selfacetbl);
    selfacetbl = crossIndexArray(selfacetbl, newoldfacetbl, {'oldfaces'});
    
    newfaces = selfacetbl.get('newfaces');
    
end
    
    

function tabcelltbl = getTabCellTbl(cccelltbl, cclinewidth, indj0, tabparams, windingnumbers, data)

    nas     = data.nas;
    nrs     = data.nrs;
    celltbl = data.celltbl;
    G       = data.G;

    for ind = 1 : numel(windingnumbers)

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

        % We can compute tabwidths as follows (commented because not used)
        % tabwidths(ind) = l(indl);

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

end


%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
