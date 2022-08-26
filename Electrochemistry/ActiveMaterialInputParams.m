classdef ActiveMaterialInputParams < ElectronicComponentInputParams
%
% Input parameter class for :class:`ActiveMaterial <Electrochemistry.ActiveMaterial>`
% 
    properties
        
        %
        % Input parameter for the  interface :class:`InterfaceInputParams <Electrochemistry.InterfaceInputParams>`
        %
        Interface
        %
        % Input parameter for the solid diffusion model  :class:`SolidDiffusionModelInputParams <Electrochemistry.SolidDiffusionModelInputParams>`
        %
        SolidDiffusion
        
        amName % Given name for the active material
                
        
        InterDiffusionCoefficient % Interdiffusion coefficient parameter (diffusion between the particles)
        
        thermalConductivity % Intrinsic Thermal conductivity of the active component
        heatCapacity        % Intrinsic Heat capacity of the active component

        electricalConductivity % Electrical conductivity / [S m^-1]

        externalCouplingTerm % structure to describe external coupling (used in absence of current collector)

        useSimplifiedDiffusionModel % Flag : true if we use simplified diffusion model

    end

    methods

        function paramobj = ActiveMaterialInputParams(jsonstruct)
            paramobj = paramobj@ElectronicComponentInputParams(jsonstruct);

            useSimplifiedDiffusionModel = paramobj.useSimplifiedDiffusionModel;
            
            pick = @(fd) pickField(jsonstruct, fd);

            paramobj.Interface = InterfaceInputParams(pick('Interface'));
           
            if useSimplifiedDiffusionModel
                paramobj.SolidDiffusion = SimplifiedSolidDiffusionModelInputParams(pick('SolidDiffusion'));                
            else
                paramobj.SolidDiffusion = FullSolidDiffusionModelInputParams(pick('SolidDiffusion'));
            end
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
