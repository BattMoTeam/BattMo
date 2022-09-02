classdef SolidElectrodeInterfaceInputParams < InputParams 

    properties
        
        molecularWeight % SEI molecular weight [kg/mol]
        density         % SEI densisity [kg/m^3]
        D               % SEI diffusion coefficient [m^2/s]
        N               % Number of discretization intervals in the sei layer model [-]
        np              % Number of computational grid cells (typically set by parent model :class:`ActiveMaterial <Electrochemistry.ActiveMaterial>`)
    end
    
    methods
        
        function paramobj = SolidElectrodeInterfaceInputParams(jsonstruct)

            paramobj = paramobj@InputParams(jsonstruct);
            
        end
        
    end
    
    
end



%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
