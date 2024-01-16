classdef TwoPointFiniteVolumeGeometry

    properties (SetAccess = immutable)

        % cell data with fields
        % - cells.centroids
        % - cells.volumes 
        cells 

        % face data with fields
        % - faces.centroids
        % - faces.areas
        % - faces.normals
        faces
        
        % node data with fields
        % - nodes.coords
        nodes

        % half-transmissibilities on each pair of cell-face
        hT

        % Only for 1d model
        faceArea
        
    end

    methods

        function tPFVgeometry = TwoPointFiniteVolumeGeometry(tPFVgeometryInput)
        % create handle

            tPFVgeometry.cells = tPFVgeometryInput.cells;
            tPFVgeometry.faces = tPFVgeometryInput.faces;
            tPFVgeometry.nodes = tPFVgeometryInput.nodes;
            tPFVgeometry.hT    = tPFVgeometryInput.hT;

            if isfield(tPFVgeometryInput, 'faceArea') | isprop(tPFVgeometryInput, 'faceArea')
                tPFVgeometry.faceArea = tPFVgeometryInput.faceArea;
            end
            
        end

    end
    
end

