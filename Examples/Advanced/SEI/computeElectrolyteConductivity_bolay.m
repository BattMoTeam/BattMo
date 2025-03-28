function conductivity = computeElectrolyteConductivity_bolay(c, T)

	ce = 1e-3 * c ;
            
    
    % Ionic conductivity, [S m^-1]
 
	conductivity = (3.4 * ce - 4.7 * ce.^1.5 + 2 * ce.^2) .* (1 + 0.2 * ce.^4).^-1;
    
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
