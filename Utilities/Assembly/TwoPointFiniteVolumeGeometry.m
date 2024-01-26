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

        % Transmissibility for internal faces (with uniform "diffusion" coefficient equation to one)
        T
        
        % Only for 1d model
        faceArea
        
    end

    methods

        function tPFVgeometry = TwoPointFiniteVolumeGeometry(tPFVgeometryInput)

            tPFVgeometry.cells = tPFVgeometryInput.cells;
            tPFVgeometry.faces = tPFVgeometryInput.faces;
            tPFVgeometry.nodes = tPFVgeometryInput.nodes;
            tPFVgeometry.hT    = tPFVgeometryInput.hT;
            tPFVgeometry.T     = tPFVgeometryInput.T ;
            
            if isfield(tPFVgeometryInput, 'faceArea') | isprop(tPFVgeometryInput, 'faceArea')
                tPFVgeometry.faceArea = tPFVgeometryInput.faceArea;
            end
            
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
