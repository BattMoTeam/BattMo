classdef SubGrid < GenericGrid

% API functions, see GenericGrid

    properties

        % Parent grid
        % The parent grid contains the updated TwoPointFiniteVolumeGeometry values (handle)
        parentGrid % Instance of Grid class

        % Mapping structure from sub-grid to parent grid with fields (see function removeCells) and inverse
        % mappings. The inverse mappings are set up in the constructor (they do not have to be given there in the
        % input). For the moment, we only need invcellmap. It is used to setup the coupling terms.
        % - mappings.cellmap
        % - mappings.facemap
        % - mappings.cellfacemap (mapping for cell-face indexing)
        % - mappings.intfacemap
        % - mappings.nodemap
        % - mappings.invcellmap
        mappings
        
    end

    methods

        function subG = SubGrid(G, parentGrid, mappings)
        % We initiate the subgrid using a MRST grid with the parentGrids and the mappings

            topology.cells.facePos   = G.cells.facePos;
            topology.cells.faces     = G.cells.faces;
            topology.cells.num       = G.cells.num;
            topology.faces.nodePos   = G.faces.nodePos;
            topology.faces.nodes     = G.faces.nodes;
            topology.faces.num       = G.faces.num;
            topology.faces.neighbors = G.faces.neighbors;
            topology.nodes.num       = G.nodes.num;
            topology.griddim         = G.griddim;

            subG.topology   = topology;
            subG.parentGrid = parentGrid;
            subG.mappings   = mappings;

            subG = subG.setupExtraMappings();
            subG = subG.setupHelpers();

        end


        function subG = setupExtraMappings(subG)

            pnc       = subG.parentGrid.topology.cells.num;
            nc        = subG.topology.cells.num;
            cellmap   = subG.mappings.cellmap;
            facemap   = subG.mappings.facemap;
            pintfaces = subG.parentGrid.helpers.intfaces;

            intfaces = all(subG.topology.faces.neighbors> 0, 2);
            intfaces = find(intfaces);
            
            invcellmap = zeros(pnc, 1);
            invcellmap(cellmap) = (1 : nc)';

            pnc     = subG.parentGrid.topology.cells.num;
            nc      = subG.topology.cells.num;
            cellmap = subG.mappings.cellmap;

            invcellmap = zeros(pnc, 1);
            invcellmap(cellmap) = (1 : nc)';

            intfacefacetbl.faces    = intfaces;
            intfacefacetbl.intfaces = (1 : numel(intfaces))';
            intfacefacetbl = IndexArray(intfacefacetbl);
            
            facepfacetbl.pfaces = facemap;
            facepfacetbl.faces = (1 : numel(facemap))';
            facepfacetbl = IndexArray(facepfacetbl);

            pintfacepfacetbl.pfaces = pintfaces;
            pintfacepfacetbl.pintfaces = (1 : numel(pintfaces))';
            pintfacepfacetbl = IndexArray(pintfacepfacetbl);
            
            intfacefacepfacetbl = crossIndexArray(intfacefacetbl, facepfacetbl, {'faces'});
            intfacefacepfacepintfacetbl = crossIndexArray(intfacefacepfacetbl, pintfacepfacetbl, {'pfaces'});

            intfacepintfacetbl = sortIndexArray(intfacefacepfacepintfacetbl, {'intfaces', 'pintfaces'});

            intfacemap = intfacepintfacetbl.get('pintfaces');

            cellpcelltbl.pcells = cellmap;
            cellpcelltbl.cells = (1 : numel(cellmap));
            cellpcelltbl = IndexArray(cellpcelltbl);

            tbls = setupTables(subG.topology);
            cellfacetbl = tbls.cellfacetbl;
            cellfacetbl = cellfacetbl.addInd('cellfaces', (1 : cellfacetbl.num)');
            
            tbls = setupTables(subG.parentGrid.topology);
            pcellpfacetbl = tbls.cellfacetbl;
            pcellpfacetbl = replacefield(pcellpfacetbl, {{'faces', 'pfaces'}, {'cells', 'pcells'}});
            pcellpfacetbl = pcellpfacetbl.addInd('pcellfaces', (1 : pcellpfacetbl.num)');

            % We take cross products to obtain the pairing between cellfaces and pcellfaces.
            tbl = crossIndexArray(cellfacetbl, facepfacetbl, {'faces'});
            tbl = crossIndexArray(tbl, cellpcelltbl, {'cells'});
            tbl = crossIndexArray(tbl, pcellpfacetbl, {'pcells', 'pfaces'});
            
            tbl = sortIndexArray(tbl, {'cellfaces', 'pcellfaces'});

            cellfacemap = tbl.get('pcellfaces');

            subG.mappings.intfacemap  = intfacemap;
            subG.mappings.cellfacemap = cellfacemap;
            subG.mappings.invcellmap  = invcellmap;
            
        end

        function subG = setupHelpers(subG)

            ptbls = setupTables(subG.parentGrid.topology);
            pcellpfacetbl = ptbls.cellfacetbl;

            tbls  = setupTables(subG.topology, 'includetbls', {'cellintfacetbl', 'extfacetbl'});
            celltbl        = tbls.celltbl;
            facetbl        = tbls.facetbl;
            extfacetbl     = tbls.extfacetbl;
            intfacetbl     = tbls.intfacetbl;
            cellintfacetbl = tbls.cellintfacetbl;
            cellfacetbl    = tbls.cellfacetbl;
            
            sgn = ones(cellintfacetbl.num, 1);
            f = cellintfacetbl.get('faces');
            c = cellintfacetbl.get('cells');
            sgn(subG.topology.faces.neighbors(f, 2) == c) = -1;

            prod = TensorProd();
            prod.tbl1 = cellintfacetbl;
            prod.tbl2 = intfacetbl;
            prod.tbl3 = celltbl;
            prod.reducefds = {'faces'};
            prod = prod.setup;

            divM = SparseTensor();
            divM = divM.setFromTensorProd(sgn, prod);
            divM = divM.getMatrix();

            diffop.div  = divM;
            diffop.grad = -divM';

            cellpcelltbl = celltbl.addInd('pcells', subG.mappings.cellmap);
            facepfacetbl = facetbl.addInd('pfaces', subG.mappings.facemap);

            pcellpfacetbl = replacefield(pcellpfacetbl, {{'cells', 'pcells'}, {'faces', 'pfaces'}});

            cellpcellpfacetbl     = crossIndexArray(pcellpfacetbl, cellpcelltbl, {'pcells'});
            cellpcellfacepfacetbl = crossIndexArray(cellpcellpfacetbl, facepfacetbl, {'pfaces'});

            % We need the following mapping as the creation of cellpcellfacepfacetbl do not guarantee that index
            % ordering is preserved.
            map1 = TensorMap();
            map1.fromTbl = pcellpfacetbl;
            map1.toTbl = cellpcellfacepfacetbl;
            map1.mergefds = {'pcells', 'pfaces'};
            map1 = map1.setup();

            P1 = map1.getMatrix();

            map2 = TensorMap();
            map2.fromTbl = cellpcellfacepfacetbl;
            map2.toTbl = cellintfacetbl;
            map2.mergefds = {'cells', 'faces'};
            map2 = map2.setup();

            P2 = map2.getMatrix();

            P = P2*P1;

            map = TensorMap();
            map.fromTbl = celltbl;
            map.toTbl = cellintfacetbl;
            map.mergefds = {'cells'};
            map = map.setup();

            D = map.getMatrix();

            map = TensorMap();
            map.fromTbl = cellintfacetbl;
            map.toTbl = intfacetbl;
            map.mergefds = {'faces'};
            map = map.setup();

            S = map.getMatrix();

            trans = struct('P', P, ...
                           'D', D, ...
                           'S', S);

            extfacepfacetbl          = crossIndexArray(extfacetbl, facepfacetbl, {'faces'});
            extfacepcellpfacetbl     = crossIndexArray(extfacepfacetbl, pcellpfacetbl, {'pfaces'});
            cellextfacepcellpfacetbl = crossIndexArray(extfacepcellpfacetbl, cellpcelltbl, {'pcells'});

            extfaces.faces = cellextfacepcellpfacetbl.get('faces');
            extfaces.cells = cellextfacepcellpfacetbl.get('cells');

            sgn = ones(cellextfacepcellpfacetbl.num, 1);
            f = extfaces.faces;
            c = extfaces.cells;
            sgn(subG.topology.faces.neighbors(f, 2) == c) = -1;

            extfaces.sgn = sgn;

            map = TensorMap();
            map.fromTbl  = cellfacetbl;
            map.toTbl    = cellextfacepcellpfacetbl;
            map.mergefds = {'cells', 'faces'};

            extfaces.cellfacemap = map.getDispatchInd();

            faceextfacemap = zeros(facetbl.num, 1);
            faceextfacemap(cellextfacepcellpfacetbl.get('faces')) = (1 : cellextfacepcellpfacetbl.num)';

            intfaces = intfacetbl.get('faces');

            subG.helpers = struct('diffop'        , diffop        , ...
                                  'trans'         , trans         , ...
                                  'extfaces'      , extfaces      , ...
                                  'faceextfacemap', faceextfacemap, ...
                                  'intfaces'      , intfaces);

        end

        function subG = setupCellFluxOperators(subG)

            G = subG.mrstFormat();
            subG.cellFluxOperators = getCellFluxOperatorsAll(G);

        end

        function G = mrstFormat(subG)

            G  = subG.topology;
            m  = subG.mappings;
            tg = subG.parentGrid.tPFVgeometry();

            d = G.griddim;

            G.cells.centroids = reshape(tg.cells.centroids, d, [])';
            G.cells.centroids = G.cells.centroids(m.cellmap, :);
            G.cells.volumes   = tg.cells.volumes(m.cellmap);

            G.faces.centroids = reshape(tg.faces.centroids, d, [])';
            G.faces.centroids = G.faces.centroids(m.facemap, :);
            G.faces.normals   = reshape(tg.faces.normals, d, [])';
            G.faces.normals   = G.faces.normals(m.facemap, :);
            G.faces.areas     = tg.faces.areas(m.facemap);

            G.nodes.coords = reshape(tg.nodes.coords, d, [])';
            G.nodes.coords = G.nodes.coords(m.nodemap, :);

            G.type = 'generic';

        end


        function vols = getVolumes(grid)

            tpfvGeometry = grid.getTPFVgeometry();
            vols = tpfvGeometry.cells.volumes;
            cmap = grid.mappings.cellmap;

            vols = vols(cmap);

        end

        function tpfvGeometry = getTPFVgeometry(grid)
            tpfvGeometry = grid.parentGrid.tPFVgeometry;
        end
                
        function areas = getFaceAreas(grid)

            tpfvGeometry = grid.getTPFVgeometry();
            areas = tpfvGeometry.faces.areas;
            fmap  = grid.mappings.facemap;

            areas = areas(fmap);

        end

        function [bchT, bccells, bcsgn] = getBcTrans(grid, bcfaces)
        % Returns half transmissibilities weighted with values of c and cell indexing for the given boundary faces
        % (bcfaces is given using subgrid indexing)
            
            tpfvGeometry = grid.getTPFVgeometry();
            hT           = tpfvGeometry.hT;
            cfmap        = grid.mappings.cellfacemap;
            exf          = grid.helpers.extfaces;

            extfaceind = grid.helpers.faceextfacemap(bcfaces);

            bccells = exf.cells(extfaceind);
            bcsgn   = exf.sgn(extfaceind);

            hTind = exf.cellfacemap(extfaceind);
            hTind = cfmap(hTind);
            
            bchT = hT(hTind);

        end

        function T = getTrans(grid)

            tpfvGeometry = grid.getTPFVgeometry();
            T            = tpfvGeometry.T;
            
            intfmap = grid.mappings.intfacemap;

            T = T(intfmap);

        end
        
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
