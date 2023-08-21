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

    end

    methods

        function grid = Grid(G, varargin)
        % We initialize the grid using a MRST grid structure

            opt = struct('faceArea', 1);
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
                grid.faceArea = opt.faceArea;
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
        
    end
    
end
