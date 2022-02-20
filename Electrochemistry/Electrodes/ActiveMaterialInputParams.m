classdef ActiveMaterialInputParams < InputParams 
%
% Input parameter class for :class:`ActiveMaterial <Electrochemistry.Electrodes.ActiveMaterial>`
%    
    properties
        
        G  % Grid

        name % Given name
        
        theta0   % Lithiation value at 0% SOC [-]
        theta100 % Lithiation value at 100% SOC[-]

        %
        % struct with fields
        %
        %  * cmax            % [mol m^-3]
        %  * D0              % [m^2 s^-1]
        %  * EaD             % [J mol^-1]
        Li 
        
        k0                      % [m^2.5 mol^-0.5 s^-1]
        Eak                     % Activation energy [J mol^-1]
        rp                      % Particle radius [m]
        volumetricSurfaceArea   % Volumetric surface area [m2 m^-3]
        volumeFraction          % Volume fraction of the active material
        density                 % Density of the active material [kg m^-3]
        n                       % number of electron transfer

        updateOCPFunc % Function to update OCP value (matlab function handler)
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
