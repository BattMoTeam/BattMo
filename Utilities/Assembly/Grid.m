classdef Grid

    properties

        % Description of topology using MRST grid structures: The fieds are
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

        % Instance of TwoPointFiniteVolumeGeometry (handle)
        % Using topology and the values of the node coordinates (and faceArea in 1D), we can compute all the other
        % properties of tPFVgeometry
        tPFVgeometry

        % Helper structures that are used in the computation of the gradient, divergence, harmonic face averages, ...
        % - helpers.diffop.grad                 (sparse matrix used in getGradient)
        % - helpers.diffop.div                  (sparse matrix used in getDiv)
        % - helpers.trans.D                     (sparse matrix used in getHarmFace method)
        % - helpers.trans.P                     (sparse matrix used in getHarmFace method)
        % - helpers.trans.S                     (sparse matrix used in getHarmFace method)
        % - helpers.extfaces.faces              (index of the external faces, sub-grid indexing)
        % - helpers.extfaces.cells              (index of the corresponding cells, sub-grid indexing)
        % - helpers.extfaces.halfTransParentInd (index of the corresponding half-transmissibility values in parent grid indexing)
        % - helpers.faceextfacemap              (mapping from face to extface, sub-grid indexing)
        helpers

        % Operators to compute norm of the flux velocity at the cell centers, see getCellFluxNorm
        % - cellFluxOperators.P
        % - cellFluxOperators.S
        cellFluxOperators
        
    end

    methods

        function grid = Grid(G, varargin)
        % We initialize the grid using a MRST grid structure

            opt = struct('faceArea', []);
            opt = merge_options(opt, varargin{:});
            
            topology.cells.facePos   = G.cells.facePos;
            topology.cells.faces     = G.cells.faces;
            topology.cells.num       = G.cells.num;
            topology.faces.nodePos   = G.faces.nodePos;
            topology.faces.nodes     = G.faces.nodes;
            topology.faces.num       = G.faces.num;
            topology.faces.neighbors = G.faces.neighbors;
            topology.nodes.num       = G.nodes.num;
            topology.griddim         = G.griddim;
            
            grid.topology   = topology;

            tf = TwoPointFiniteVolumeGeometry();

            if topology.griddim == 1
                if isempty(opt.faceArea)
                    tf.faceArea = 1;
                else
                    tf.faceArea = opt.faceArea;
                end
            end

            % We flatten node coordinates structure
            tf.nodes.coords = reshape(G.nodes.coords', [], 1);
            
            grid.tPFVgeometry = tf;
            
            grid = grid.setupHelpers();
            
        end

        function updateTPFgeometry(grid)

            G  = grid.topology;
            d  = grid.topology.griddim;
            tf = grid.tPFVgeometry;
            
            if  G.griddim == 1
                farea = tf.faceArea;
            else
                farea = 1;
            end
            
            G.nodes.coords = reshape(tf.nodes.coords, d, [])';

            G.type = 'generic';
            G = computeGeometry(G);
            
            rock.poro = ones(G.cells.num, 1);
            rock.perm = ones(G.cells.num, 1);
            hT = computeTrans(G, rock);

            cells.centroids = reshape(G.cells.centroids', [], 1);
            cells.volumes   = farea*G.cells.volumes;
            
            faces.centroids = reshape(G.faces.centroids', [], 1);
            faces.normals   = farea*reshape(G.faces.normals', [], 1);
            faces.areas     = farea*G.faces.areas;

            grid.tPFVgeometry.cells = cells; 
            grid.tPFVgeometry.faces = faces;
            grid.tPFVgeometry.hT    = farea*hT;

        end

        function grid = setupHelpers(grid)
        % setup helpers

            tp = grid.topology;
            
            tbls = setupTables(tp, 'includetbls', {'intfacetbl', 'extfacetbl'});
            intfacetbl     = tbls.intfacetbl;
            cellintfacetbl = tbls.cellintfacetbl;
            cellfacetbl    = tbls.cellfacetbl;
            celltbl        = tbls.celltbl;
            facetbl        = tbls.facetbl;
            extfacetbl     = tbls.extfacetbl;

            sgn = ones(cellintfacetbl.num, 1);
            f = cellintfacetbl.get('faces');
            c = cellintfacetbl.get('cells');
            sgn(grid.topology.faces.neighbors(f, 1) == c) = -1;

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
            
            map = TensorMap();
            map.fromTbl  = cellfacetbl;
            map.toTbl    = cellintfacetbl;
            map.mergefds = {'cells', 'faces'};
            map = map.setup();
            
            P = map.getMatrix();

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
            
            cellextfacetbl = crossIndexArray(cellfacetbl, extfacetbl, {'faces'});

            extfaces.faces = cellextfacetbl.get('faces');
            extfaces.cells = cellextfacetbl.get('cells');

            map = TensorMap();
            map.fromTbl  = cellfacetbl;
            map.toTbl    = cellextfacetbl;
            map.mergefds = {'cells', 'faces'};
            
            extfaces.halfTransParentInd = map.getDispatchInd();

            sgn = ones(cellextfacetbl.num, 1);
            f = extfaces.faces;
            c = extfaces.cells;
            sgn(grid.topology.faces.neighbors(f, 1) == c) = -1;

            extfaces.sgn = sgn;

            faceextfacemap = zeros(facetbl.num, 1);
            faceextfacemap(extfacetbl.get('faces')) = (1 : extfacetbl.num)';

            grid.helpers = struct('diffop'        , diffop        , ...
                                  'trans'         , trans         , ...
                                  'extfaces'      , extfaces      , ...
                                  'faceextfacemap', faceextfacemap);
        end

        function G = getMRSTgrid(grid)

            G  = grid.topology;
            tg = grid.tPFVgeometry;

            d = G.griddim;

            G.type = 'generic';

            G.cells.centroids = reshape(tg.cells.centroids, d, [])';
            G.cells.volumes   = tg.cells.volumes;
            G.faces.centroids = reshape(tg.faces.centroids, d, [])';
            G.faces.normals   = reshape(tg.faces.normals, d, [])';
            G.faces.areas     = tg.faces.areas;
            G.nodes.coords    = reshape(tg.nodes.coords, d, [])';
            
        end

        function grid = setupCellFluxOperators(grid)

            grid.updateTPFgeometry();
            G = getMRSTgrid(grid);
            grid.cellFluxOperators = getCellFluxOperatorsAll(G);
            
        end
        
        function vols = getVolumes(grid)
            
            vols = grid.tPFVgeometry.cells.volumes;
            
        end

        function areas = getFaceAreas(grid)
            
            areas = grid.tPFVgeometry.faces.areas;
            
        end
        
        function u = getHarmFace(grid, c)
        % Returns fluxes for each internal faces for the cell-valued vector c

            op = grid.helpers.trans;
            hT = grid.tPFVgeometry.hT;
            
            u = 1 ./ (op.S * ( 1 ./ (op.D*c .* op.P*hT)));
            
        end
        
        function [bchT, bccells, bcsgn] = getBcHarmFlux(grid, u, bcfaces)
        % Returns half transmissibilities and cell indexing for the given boundary faces

            hT   = grid.tPFVgeometry.hT;
            exf  = grid.helpers.extfaces;

            extfaceind = grid.helpers.faceextfacemap(bcfaces);
            
            bccells = exf.cells(extfaceind);
            bcsgn   = exf.sgn(extfaceind);
            bchT    = hT(exf.halfTransParentInd(extfaceind));
            
        end

        function jsq = getCellFluxNorm(grid, u)
    
            P = grid.cellFluxOperators.P;
            S = grid.cellFluxOperators.S;
    
            j = P*u;
            jsq = j.^2;
            jsq = S*jsq;
            
        end

        function nc = getNumberOfCells(grid)

            nc = grid.topology.cells.num;
            
        end        
    end
    
end
