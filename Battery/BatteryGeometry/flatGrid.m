function output = flatGrid(params)
    
    
    %% Parameters 
    % 
    % params structure with following fields
    % - nwindings : number of windings in the spiral
    % - depth     : depth of domain
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
    depth     = params.depth;
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

    h = linspace(0, depth, nas + 1);

    nperlayer = sum(nrs);

    G = tensorGrid(h, w);
    G.faces = rmfield(G.faces, 'tag');
    G.type = {'flatGrid'};
    G = computeGeometry(G, 'findNeighbors', true);

    n = numel(h) - 1;
    m = numel(w) - 1;

    ncomp = numel(widths);
    comptag = rldecode((1 : ncomp)', nrs);
    comptag = repmat(comptag, [nwindings, 1]);

    comptagtbl.tag = comptag;
    comptagtbl.indj = (1 : (sum(nrs)*nwindings))';
    comptagtbl = IndexArray(comptagtbl);

    celltbl.cells = (1 : G.cells.num)';
    celltbl.indi = repmat((1 : nas)', [sum(nrs)*nwindings, 1]);
    celltbl.indj = rldecode((1 : sum(nrs)*nwindings)', nas*ones(sum(nrs)*nwindings, 1));
    celltbl = IndexArray(celltbl);

    celltagtbl = crossIndexArray(celltbl, comptagtbl, {'indj'});
    celltagtbl = sortIndexArray(celltagtbl, {'cells', 'tag'});

    tag = celltagtbl.get('tag');
    
    %% Extrude battery in z-direction
    
    zwidths = (L/nL)*ones(nL, 1);
    if nL == 1
        G = makeLayeredGrid(G, 1);
        k  = G.nodes.coords(:, 3) > 0;
        G.nodes.coords(k, 3) = zwidths;    
    else
        G = makeLayeredGrid(G, zwidths);
    end
    G = computeGeometry(G);
    
    tag = repmat(tag, [nL, 1]);
    
    tbls = setupSimpleTables(G);
    cellfacetbl = tbls.cellfacetbl;
    
    clear extfacetbl
    extfacetbl.faces = find(any(G.faces.neighbors == 0, 2));
    extfacetbl = IndexArray(extfacetbl);
    extcellfacetbl = crossIndexArray(extfacetbl, cellfacetbl, {'faces'});

    %% We set up Cartesian indexing for the cells

    [indi, indj, indk] = ind2sub([n, m, nL], (1 : G.cells.num)');
    
    clear celltbl
    celltbl.cells = (1 : G.cells.num)';
    celltbl.indi = indi;
    celltbl.indj = indj;
    celltbl.indk = indk;
    celltbl = IndexArray(celltbl);

    %% We add vertical (1) and horizontal (2) direction index for the faces (see makeLayeredGrid for the setup)
    
    nf = G.faces.num;
    clear facetbl
    facetbl.faces = (1 : nf)';
    dir = 2*ones(nf, 1);
    dir(1 : (nL + 1)*n*m) = 1;
    facetbl.dir = dir;
    facetbl = IndexArray(facetbl);

    
    %% We extract the faces at the exterior for thermal exchange, using Cartesian indexing
    
    scelltbl.indi = (1: n)';
    scelltbl.indj = m*ones(n, 1);
    scelltbl      = IndexArray(scelltbl);
    
    scelltbl = crossIndexArray(celltbl, scelltbl, {'indi', 'indj'});
    scelltbl = projIndexArray(scelltbl, 'cells');
    extscellfacetbl = crossIndexArray(scelltbl, extcellfacetbl, {'cells'});
    sfacetbl = projIndexArray(extscellfacetbl, {'faces'});
    
    sfacetbl = sfacetbl.addInd('dir', 2*ones(sfacetbl.num, 1));
    sfacetbl = crossIndexArray(sfacetbl, facetbl, {'faces', 'dir'});
    
    sfaces = sfacetbl.get('faces');

    clear scelltbl
    nnrs = sum(nrs);
    scelltbl.indk = [1; nL];
    scelltbl = IndexArray(scelltbl);
    
    scelltbl = crossIndexArray(celltbl, scelltbl, {'indk'});
    scelltbl = projIndexArray(scelltbl, 'cells');
    extscellfacetbl = crossIndexArray(scelltbl, extcellfacetbl, {'cells'});
    sfacetbl = projIndexArray(extscellfacetbl, {'faces'});
    
    sfacetbl = sfacetbl.addInd('dir', 1*ones(sfacetbl.num, 1));
    sfacetbl = crossIndexArray(sfacetbl, facetbl, {'faces', 'dir'});

    sfaces = [sfaces; sfacetbl.get('faces')];
    
    thermalExchangeFaces = sfaces;

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
