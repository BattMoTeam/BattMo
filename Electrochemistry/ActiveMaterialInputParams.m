classdef ActiveMaterialInputParams < ElectroChemicalComponentInputParams
%
% Input parameter class for :class:`ActiveMaterial <Electrochemistry.ActiveMaterial>`
% 
    properties
        
        %
        % Input parameter for the active material  :class:`InterfaceInputParams <Electrochemistry.Electrodes.InterfaceInputParams>`
        %
        Interface
        
        amName % Given name for the active material
                
        
        InterDiffusionCoefficient % Interdiffusion coefficient parameter (diffusion between the particles)
        
        thermalConductivity % Intrinsic Thermal conductivity of the active component
        heatCapacity       % Intrinsic Heat capacity of the active component

        electricalConductivity % Electrical conductivity / [S m^-1]

    end

    methods

        function paramobj = ActiveMaterialInputParams(jsonstruct)
            paramobj = paramobj@ElectroChemicalComponentInputParams(jsonstruct);

            pick = @(fd) pickField(jsonstruct, fd);
            paramobj.Interface = InterfaceInputParams(pick('Interface'));
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
