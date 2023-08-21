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
        
        % Value of the node coordinates
        nodecoords

        % Value of the face area - only for 1D model
        faceArea

        % Instance of TwoPointFiniteVolumeGeometry (handle)
        % Using topology and nodecoords all the properties of tPFVgeometry can be updated
        tPFVgeometry

        % Helper structures that are used to extract the half-transmissibilities from the parent structures and assemble the
        % fluxes, with fields
        % - helpers.trans.D                     (sparse matrix used in getFlux method)
        % - helpers.trans.P                     (sparse matrix used in getFlux method)
        % - helpers.trans.S                     (sparse matrix used in getFlux method)
        % - helpers.extfaces.faces              (index of the external faces, sub-grid indexing)
        % - helpers.extfaces.cells              (index of the corresponding cells, sub-grid indexing)
        % - helpers.extfaces.halfTransParentInd (index of the corresponding half-transmissibility values in parent grid indexing)
        % - helpers.faceextfacemap              (mapping from face to extface, sub-grid indexing)
        helpers
        
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

            if topology.griddim == 1
                if isempty(opt.faceArea)
                    grid.faceArea = 1;
                else
                    grid.faceArea = opt.faceArea;
                end
            end

            % We flatten the coordinate structures
            grid.nodecoords = reshape(G.nodes.coords', [], 1);
            
            grid.tPFVgeometry = TwoPointFiniteVolumeGeometry();
            
        end

        function updateTPFgeometry(grid)

            G          = grid.topology;
            nodecoords = grid.nodecoords;
            d          = grid.topology.griddim;
            farea      = grid.faceArea;
            
            G.nodes.coords = reshape(nodecoords, d, [])';

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

            nodes.coords = nodecoords;

            grid.tPFVgeometry.cells = cells; 
            grid.tPFVgeometry.faces = faces;
            grid.tPFVgeometry.nodes = nodes;
            grid.tPFVgeometry.hT    = farea*hT;

            % setup helpers

            tbls = setupTables(G, 'includetbls', {'intfacetbl', 'extfacetbl'});
            intfacetbl     = tbls.intfacetbl;
            cellintfacetbl = tbls.cellintfacetbl;
            cellfacetbl    = tbls.cellfacetbl;
            celltbl        = tbls.celltbl;
            extfacetbl     = tbls.extfacetbl;
            
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

            map = TensorMap();
            map.fromTbl  = cellfacetbl;
            map.toTbl    = cellextfacetbl;
            map.mergefds = {'cells', 'faces'};
            
            extfaces.faces              = cellextfacetbl.get('faces');
            extfaces.cells              = cellextfacetbl.get('cells');
            extfaces.halfTransParentInd = map.getDispatchInd();

            faceextfacemap = zeros(facetbl.num, 1);
            faceextfacemap(extfacetbl.get('faces')) = (1 : extfacetbl.num)';

            grid.helpers = struct('trans'         , trans   , ...
                                  'extfaces'      , extfaces, ...
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

        function vols = getVolumes(grid)
            
            vols = grid.tPFVgeometry.cells.volumes;
            
        end
        
        function u = getFlux(grid, c)
        % Returns fluxes for each internal faces for the cell-valued vector c

            op = grid.helpers.trans;
            hT = grid.tPFVgeometry.hT;
            
            u = 1 ./ (op.S * ( 1 ./ (op.D*c .* op.P*hT)));
            
        end
        
        function [bchT, bccells] = getBcFlux(grid, u, bcfaces)
        % Returns half transmissibilities and cell indexing for the given boundary faces

            hT   = grid.tPFVgeometry.hT;
            exf  = grid.helpers.extfaces;

            extfaceind = grid.helpers.faceextfacemap(bcfaces);
            
            bccells = exf.cells(extfaceind);
            bchT    = hT(exf.halfTransParentInd(extfaceind));

        end
        
    end
    
end
