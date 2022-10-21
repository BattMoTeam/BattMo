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
        
        InterDiffusionCoefficient % Interdiffusion coefficient parameter (diffusion between the particles)
        
        thermalConductivity % Intrinsic Thermal conductivity of the active component
        heatCapacity        % Intrinsic Heat capacity of the active component

        electricalConductivity % Electrical conductivity / [S m^-1]

        externalCouplingTerm % structure to describe external coupling (used in absence of current collector)

        diffusionModelType % Choose between type of diffusion model ('full' or 'simple'. The default is set to 'full')

        BruggemanCoefficient

        volumeFraction % Volume fraction of the whole material (binder and so on included)
        
        activeMaterialFraction = 1 % Volume fraction occupied only by the active material (default value is 1)

        use_particle_diffusion

    end

    methods

        function paramobj = ActiveMaterialInputParams(jsonstruct)

            paramobj = paramobj@ElectronicComponentInputParams(jsonstruct);

            if isempty(paramobj.diffusionModelType)
                % we do not use any diffusion model (use_particle_diffusion = false)
            end
            
            pick = @(fd) pickField(jsonstruct, fd);

            paramobj.Interface = InterfaceInputParams(pick('Interface'));

            diffusionModelType = paramobj.diffusionModelType;

            switch diffusionModelType
                
              case 'simple'
                
                paramobj.SolidDiffusion = SimplifiedSolidDiffusionModelInputParams(pick('SolidDiffusion'));
                
              case 'full'

                paramobj.SolidDiffusion = FullSolidDiffusionModelInputParams(pick('SolidDiffusion'));
                
              otherwise
                
                error('Unknown diffusionModelType %s', diffusionModelType);
                
            end

            paramobj = paramobj.validateInputParams();
            
        end

        function paramobj = validateInputParams(paramobj)


            diffusionModelType = paramobj.diffusionModelType;
            
            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            switch diffusionModelType
                
              case 'simple'
                
                paramobj = mergeParameters(paramobj, {'volumeFraction'}, {itf, 'volumeFraction'});
                
              case 'full'
                
                paramobj = mergeParameters(paramobj, {'volumeFraction'}, {itf, 'volumeFraction'});
                paramobj = mergeParameters(paramobj, {'volumeFraction'}, {sd, 'volumeFraction'});
                paramobj = mergeParameters(paramobj, {'activeMaterialFraction'}, {sd, 'activeMaterialFraction'});
                
                if ~isempty(paramobj.(sd).D)
                    % we impose that cmax in the solid diffusion model and the interface are consistent
                    paramobj = mergeParameters(paramobj, {sd, 'cmax'}, {itf, 'cmax'}, 'force', false);
                    paramobj = mergeParameters(paramobj, {sd, 'theta0'}, {itf, 'theta0'}, 'force', false);
                    paramobj = mergeParameters(paramobj, {sd, 'theta100'}, {itf, 'theta100'}, 'force', false);
                end

              otherwise
                error('Unknown diffusionModelType %s', diffusionModelType);
            end

            paramobj = validateInputParams@ElectronicComponentInputParams(paramobj);
            
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
