classdef ActiveMaterialInputParams < InputParams 
%
% Input class for :class:`ActiveMaterial <Electrochemistry.Electrodes.ActiveMaterial>`
%    
    properties
        
        G

        name
        
        theta0                  % at 0% SOC [-]
        theta100                % at 100% SOC[-]
        Li % struct with fields
           %  - cmax            % [mol m^-3]
           %  - D0              % [m^2 s^-1]
           %  - EaD             % [J mol^-1]
        k0                      % [m^2.5 mol^-0.5 s^-1]
        Eak                     % [J mol^-1]
        rp                      % [m]
        volumetricSurfaceArea   % [m2 m^-3]
        volumeFraction         
        density                 % [kg m^-3]
        n                       % number of electron transfer

        updateOCPFunc % Function to update OCP value
    end
    
    methods
        
        function paramobj = ActiveMaterialInputParams(jsonstruct);
            paramobj = paramobj@InputParams(jsonstruct);
        end
        
    end
    
    
end



%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BatMo

BatMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BatMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BatMo.  If not, see <http://www.gnu.org/licenses/>.
%}
