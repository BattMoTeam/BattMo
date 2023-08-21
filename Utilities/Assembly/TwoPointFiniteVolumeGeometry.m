classdef TwoPointFiniteVolumeGeometry < handle

    properties

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
    
end

