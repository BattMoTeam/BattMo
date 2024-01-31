classdef SeparatorInputParams < ComponentInputParams
%
% Input parameter class for :code:`Separator` model
%    
    properties
        
        porosity             % the ratio of the volume free space to the total volume (symbol: varepsilon)
        density              % the mass density of the material (symbol: rho)
        bruggemanCoefficient % coefficient to determine effective transport parameters in porous media (symbol: beta)
        
        thermalConductivity  % Intrinsic Thermal conductivity of the electrolyte
        specificHeatCapacity % Specific Heat capacity of the electrolyte

        % Advanced parameters
        effectiveThermalConductivity    % (account for volume fraction)
        effectiveVolumetricHeatCapacity % (account for volume fraction and density)

        % helper parameters
        use_thermal
    end
    
    methods

        function inputparams = SeparatorInputParams(jsonstruct)
            inputparams = inputparams@ComponentInputParams(jsonstruct);
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
