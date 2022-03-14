function output = radialGrid(params)
    
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

    h = linspace(0, 2*pi*rInner, nas + 1);

    nperlayer = sum(nrs);

    cartG = tensorGrid(h, w);

    n = numel(h) - 1;
    m = numel(w) - 1;

    % plotGrid(cartG)

    % We roll the domain into a spirale
    x = cartG.nodes.coords(:, 1);
    y = cartG.nodes.coords(:, 2);

    theta = x./rInner;

    cartG.nodes.coords(:, 1) = (rInner + y).*cos(theta);
    cartG.nodes.coords(:, 2) = (rInner + y).*sin(theta);
    
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
    node2tbl.indi1 = ones(m + 1, 1);
    node2tbl.indj1 = (1 : (m + 1))';
    node2tbl.indi2 = (n + 1)*ones(m + 1, 1);
    node2tbl.indj2 = (1 : (m + 1))';
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
    face2tbl.indi1 = ones(m, 1);
    face2tbl.indj1 = (1 : m)';
    face2tbl.indi2 = (n + 1)*ones(m , 1);
    face2tbl.indj2 = (1 : m)';
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
    G.type = {'radialGrid'};
    G = computeGeometry(G, 'findNeighbors', true);
    
    ncomp = numel(widths);
    comptag = rldecode((1 : ncomp)', nrs);
    comptag = repmat(comptag, [nwindings, 1]);

    comptagtbl.tag = comptag;
    comptagtbl.indj = (1 : (sum(nrs)*nwindings))';
    comptagtbl = IndexArray(comptagtbl);

    celltbl.cells = (1 : cartG.cells.num)';
    celltbl.indi = repmat((1 : nas)', [sum(nrs)*nwindings, 1]);
    celltbl.indj = rldecode((1 : sum(nrs)*nwindings)', nas*ones(sum(nrs)*nwindings, 1));
    celltbl = IndexArray(celltbl);

    celltagtbl = crossIndexArray(celltbl, comptagtbl, {'indj'});
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

    [indi, indj, indk] = ind2sub([n, m, nL], (1 : G.cells.num)');
    
    clear celltbl
    celltbl.cells = (1 : G.cells.num)';
    celltbl.indi = indi;
    celltbl.indj = indj;
    celltbl.indk = indk;
    celltbl = IndexArray(celltbl);

    % We add vertical (1) and horizontal (2) direction index for the faces (see makeLayeredGrid for the setup)
    
    nf = G.faces.num;
    clear facetbl
    facetbl.faces = (1 : nf)';
    dir = 2*ones(nf, 1);
    dir(1 : (nL + 1)*n*m) = 1;
    facetbl.dir = dir;
    facetbl = IndexArray(facetbl);
    
    scelltbl.indi = (1: n)';
    scelltbl.indj = 1*ones(n, 1);
    scelltbl = IndexArray(scelltbl);
    
    scelltbl = crossIndexArray(celltbl, scelltbl, {'indi', 'indj'});
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
