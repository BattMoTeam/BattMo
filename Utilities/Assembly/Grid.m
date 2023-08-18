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

        % Instance of TwoPointFiniteVolumeGeometry (handle)
        % Using topology and nodecoords all the properties of tPFVgeometry can be updated
        tPFVgeometry

    end


    methods

        function grid = Grid(G)
        % We initialize the grid using a MRST grid structure

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

            % We flatten the coordinate structures
            grid.nodecoords = reshape(G.nodes.coords', [], 1);
            
            grid.tPFVgeometry = TwoPointFiniteVolumeGeometry();
            
        end

        function updateTPFgeometry(grid)
            
        end

        function G = getMRSTgrid(grid)

        end
        
    end
    
end
