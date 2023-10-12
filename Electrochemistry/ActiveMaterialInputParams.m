classdef ActiveMaterialInputParams < ComponentInputParams
%
% Input parameter class for :code:`ActiveMaterial` model
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
        
        %% Standard parameters

        electronicConductivity % the electronic conductivity of the material (symbol: sigma)
        density                % the mass density of the material (symbol: rho)
        massFraction           % the ratio of the mass of the material to the total mass of the phase or mixture (symbol: gamma)
        
        thermalConductivity    % the intrinsic Thermal conductivity of the active component
        specificHeatCapacity   % Specific Heat capacity of the active component

        diffusionModelType     % diffusion model type, either 'full' or 'simplified'

        %% Advanced parameters

        standAlone % Set to true if Active Material is used as a stand-alone model (not within a battery cell, see runActiveMaterial example)

        %% Coupling parameters
        
        externalCouplingTerm % structure to describe external coupling (used in absence of current collector)
        
    end

    methods

        function paramobj = ActiveMaterialInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
            if isempty(paramobj.standAlone)
                paramobj.standAlone = false;
            end

            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            pick = @(fd) pickField(jsonstruct, fd);
            
            paramobj.(itf) = InterfaceInputParams(pick('Interface'));

            if isempty(paramobj.diffusionModelType)
                paramobj.diffusionModelType = 'full';
            end

            diffusionModelType = paramobj.diffusionModelType;
            
            switch diffusionModelType
                
              case 'simple'
                
                paramobj.(sd) = SimplifiedSolidDiffusionModelInputParams(pick(sd));
                
              case 'full'

                paramobj.(sd) = FullSolidDiffusionModelInputParams(pick(sd));
                
              otherwise
                
                error('Unknown diffusionModelType %s', diffusionModelType);
                
            end

            paramobj = paramobj.validateInputParams();
            
        end

        function paramobj = validateInputParams(paramobj)


            diffusionModelType = paramobj.diffusionModelType;
            
            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            paramobj = mergeParameters(paramobj, {{'density'}, {itf, 'density'}});

            if paramobj.standAlone
                % only one particle in the stand-alone model
                paramobj.(sd).np = 1;
            end
            
            switch diffusionModelType
                
              case 'simple'
                
                paramobj = mergeParameters(paramobj, {{itf, 'volumetricSurfaceArea'}, {sd, 'volumetricSurfaceArea'}});
                
              case 'full'
                
                paramobj = mergeParameters(paramobj, {{itf, 'volumetricSurfaceArea'}, {sd, 'volumetricSurfaceArea'}});
                
                if ~isempty(paramobj.(sd).diffusionCoefficient)
                    % we impose that cmax in the solid diffusion model and the interface are consistent
                    paramnames = {'saturationConcentration', ...
                                  'guestStoichiometry100', ...
                                  'guestStoichiometry0'};
                    for iparam = 1 : numel(paramnames)
                        paramname = paramnames{iparam};
                        paramobj = mergeParameters(paramobj, {{sd, paramname}, {itf, paramname}}, 'force', false);
                    end
                    
                end

              case 'interParticleOnly'

                paramobj = mergeParameters(paramobj, {{'volumeFraction'}, {itf, 'volumeFraction'}});
                paramobj.SolidDiffusion = [];
                
              otherwise
                
                error('Unknown diffusionModelType %s', diffusionModelType);
                
            end

            paramobj = validateInputParams@ComponentInputParams(paramobj);
            
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
