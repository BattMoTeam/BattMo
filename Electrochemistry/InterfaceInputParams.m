classdef InterfaceInputParams < InputParams 
%
% Input parameter class for :class:`ActiveMaterial <Electrochemistry.Electrodes.ActiveMaterial>`
%    
    properties
        
        G  % Grid
        
        theta0   % Lithiation value at 0% SOC [-]
        theta100 % Lithiation value at 100% SOC [-]

        cmax                    % Maximum concentration [mol m^-3]
        k0                      % Reference rate constant  [m^2.5 mol^-0.5 s^-1]
        Eak                     % Activation energy [J mol^-1]
        volumetricSurfaceArea   % Volumetric surface area [m2 m^-3]
        volumeFraction          % Volume fraction of the active material
        density                 % Density of the active material [kg m^-3]
        n                       % number of electron transfer

        OCP % Function to update OCP value given as a struct with fields
            % OCP.type = "function";
            % OCP.functionname :  matlab function name (should be available in path)
            % OCP.argumentlist = ["cElectrode", "T", "cmax"]

        j0  % Function to update j0 : if empty, the default expression using k0 is used, see method Interface.updateReactionRateCoefficient
            % j0.type = {"function", "constant"} % if "constant" is selected, the value of k0 is used to compute reaction rate
            % j0.functionname :  matlab function name (should be available in path)
            % j0.argumentlist = ["cElectrodeSurface", "cmax"]

        alpha = 0.5 % coefficient in Butler-Volmer coefficient (default is 0.5)
    end
    
    methods
        
        function paramobj = InterfaceInputParams(jsonstruct)

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
