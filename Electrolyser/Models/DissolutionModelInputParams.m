classdef DissolutionModelInputParams < ComponentInputParams

    properties

        E0                     % Standard electrical potential for the dissolution reaction [V]
        c0                     % Reference concentration [mol/m^3]
        j0                     % Reference exchange current density for the dissolution reaction [A/m^2]
        MW                     % Molar mass
        rho                    % Density
        volumeFraction0        % Initial volume fraction
        volumetricSurfaceArea0 % Initial volumetric surface area
        
    end

    methods
        
        function paramobj = DissolutionModelInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
end



%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
