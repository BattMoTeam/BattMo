classdef SubGrid

    properties

        % Parent grid
        % The parent grid contains the updated TwoPointFiniteVolumeGeometry values (handle)
        parentGrid % Instance of Grid class

        % Description of the topology of the sub-grid (local indexing) using MRST grid structures: The fieds are
        % - topology.cells.facePos
        % - topology.cells.faces
        % - topology.cells.num
        % - topology.faces.nodePos
        % - topology.faces.nodes
        % - topology.faces.num
        % - topology.faces.neighbors
        % - topology.nodes.num
        % - topology.griddim
        topology

        % Mapping structure from sub-grid to parent grid with fields (see function removeCells) and inverse
        % mappings. The inverse mappings are set up in the constructor (they do not have to be given there in the
        % input). For the moment, we only need invcellmap. It is used to setup the coupling terms.
        % - mappings.cellmap
        % - mappings.facemap
        % - mappings.nodemap
        % - mappings.invcellmap
        mappings

        % Helper structures that are used to extract the half-transmissibilities from the parent structures and assemble the
        % fluxes, with fields
        % - helpers.diffop.grad                 (sparse matrix used in getGradient)
        % - helpers.diffop.div                  (sparse matrix used in getDiv)
        % - helpers.trans.D                     (sparse matrix used in getTransHarmFace method)
        % - helpers.trans.P                     (sparse matrix used in getTransHarmFace method)
        % - helpers.trans.S                     (sparse matrix used in getTransHarmFace method)
        % - helpers.extfaces.faces              (index of the external faces, sub-grid indexing)
        % - helpers.extfaces.cells              (index of the corresponding cells, sub-grid indexing)
        % - helpers.extfaces.sgn                (sign of the corresponding cell-face pair)
        % - helpers.extfaces.halfTransParentInd (index of the corresponding half-transmissibility values in parent grid indexing)
        % - helpers.faceextfacemap              (mapping from face to extface, sub-grid indexing)
        % - helpers.intfaces                    (index of internal faces)
        helpers

        % Operators to compute norm of the flux velocity at the cell centers, see getCellFlux
        % - cellFluxOperators.P
        % - cellFluxOperators.S
        cellFluxOperators

        % API functions are
        % - getVolumes
        % - getFaceAreas
        % - getFaceCentroids
        % - getGrad
        % - getDiv
        % - getTrans
        % - getTransHarmFace
        % - getTransBcHarmFace
        % - getCellFluxNorm
        % - getIntFaceIndex
        % - getNumberOfCells
        % - getNumberOfFaces

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

            subG = subG.setupInverseMappings();
            subG = subG.setupHelpers();

        end


        function subG = setupInverseMappings(subG)

            pnc     = subG.parentGrid.topology.cells.num;
            nc      = subG.topology.cells.num;
            cellmap = subG.mappings.cellmap;

            invcellmap = zeros(pnc, 1);
            invcellmap(cellmap) = (1 : nc)';

            pnc     = subG.parentGrid.topology.cells.num;
            nc      = subG.topology.cells.num;
            cellmap = subG.mappings.cellmap;

            invcellmap = zeros(pnc, 1);
            invcellmap(cellmap) = (1 : nc)';

            subG.mappings.invcellmap = invcellmap;

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

            map = TensorMap();
            map.fromTbl  = pcellpfacetbl;
            map.toTbl    = cellextfacepcellpfacetbl;
            map.mergefds = {'pcells', 'pfaces'};

            extfaces.halfTransParentInd = map.getDispatchInd();

            sgn = ones(cellextfacepcellpfacetbl.num, 1);
            f = extfaces.faces;
            c = extfaces.cells;
            sgn(subG.topology.faces.neighbors(f, 2) == c) = -1;

            extfaces.sgn = sgn;

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

            G = getMRSTgrid(subG);
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
            G.faces.normals = G.faces.normals(m.facemap, :);
            G.faces.areas     = tg.faces.areas(m.facemap);

            G.nodes.coords = reshape(tg.nodes.coords, d, [])';
            G.nodes.coords = G.nodes.coords(m.nodemap, :);

            G.type = 'generic';

        end

        % API functions which should be used in simulations

        function vols = getVolumes(subG)

            pvols = subG.parentGrid.tPFVgeometry.cells.volumes;
            cmap  = subG.mappings.cellmap;

            vols = pvols(cmap);

        end

        function areas = getFaceAreas(subG)

            pareas = subG.parentGrid.tPFVgeometry.faces.areas;
            fmap  = subG.mappings.facemap;

            areas = pareas(fmap);

        end

        function c = getFaceCentroids(subG)

            d = subG.parentGrid.topology.griddim;
            c0 = reshape(subG.parentGrid.tPFVgeometry.faces.centroids, d, [])';
            fmap = subG.mappings.facemap;
            c = c0(fmap,:);

        end

        function v = getGrad(subG, c)

            v = subG.helpers.diffop.grad*c;

        end

        function u = getDiv(subG, v)

            u = subG.helpers.diffop.div*v;

        end

        function u = getTrans(subG)

            ind = subG.helpers.intfaces;
            T   = subG.parentGrid.tPFVgeometry.T;
            
            u = T(ind);
            
        end

        function u = getTransHarmFace(subG, c)
        % Returns fluxes for each internal faces for the cell-valued vector c

            t = subG.helpers.trans;
            hT = subG.parentGrid.tPFVgeometry.hT;

            u = 1 ./ (t.S * ( 1 ./ ((t.D*c) .* (t.P*hT))));

        end

        function [bchT, bccells, bcsgn] = getTransBcHarmFace(subG, c, bcfaces)
        % Returns half transmissibilities weighted with values of c and cell indexing for the given boundary faces
        % (bcfaces is given using subgrid indexing)

            hT   = subG.parentGrid.tPFVgeometry.hT;
            exf  = subG.helpers.extfaces;

            extfaceind = subG.helpers.faceextfacemap(bcfaces);

            bccells = exf.cells(extfaceind);
            bcsgn   = exf.sgn(extfaceind);
            bchT    = hT(exf.halfTransParentInd(extfaceind)).*c;

        end

        function jsq = getCellFluxNorm(subG, u)

            P = subG.cellFluxOperators.P;
            S = subG.cellFluxOperators.S;

            j = P*u;
            jsq = j.^2;
            jsq = S*jsq;

        end


        function intfaces = getIntFaceIndex(subG)

            intfaces = subG.helpers.intfaces;

        end

        function nc = getNumberOfCells(subG)

            nc = subG.topology.cells.num;

        end

        function nf = getNumberOfFaces(subG)

            nf = subG.topology.faces.num;

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
