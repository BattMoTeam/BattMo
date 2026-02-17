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

        LithiumPlating
        
        %% Standard parameters

        electronicConductivity % the electronic conductivity of the material (symbol: sigma)
        density                % the mass density of the material (symbol: rho)
        massFraction           % the ratio of the mass of the material to the total mass of the phase or mixture (symbol: gamma)
        
        thermalConductivity    % the intrinsic Thermal conductivity of the active component
        specificHeatCapacity   % Specific Heat capacity of the active component

        diffusionModelType     % diffusion model type, string equal to either
                               % - 'full'
                               % - 'simple'
                               % - 'swelling'

        %% SEI layer model choice
        SEImodel % string defining the sei model, see schema Utilities/JsonSchemas/ActiveMaterial.schema.json. Can take value
                  % - 'none' (default)
                  % - 'Bolay'
                  % - 'Safari'

        %% Advanced parameters

        isRootSimulationModel % Set to true if Active Material is used as a stand-alone model (not within a battery cell, see runActiveMaterial for an example)

        %% Coupling parameters
        
        externalCouplingTerm % structure to describe external coupling (used in absence of current collector)

        useLithiumPlating
        
    end

    methods

        function inputparams = ActiveMaterialInputParams(jsonstruct)

            sd  = 'SolidDiffusion';
            itf = 'Interface';
            lp  = 'LithiumPlating';
            
            jsonstruct = equalizeJsonStructField(jsonstruct, 'density', {itf, 'density'});

            jsonstruct = setDefaultJsonStructField(jsonstruct,  'useLithiumPlating', false);

            jsonstruct = setDefaultJsonStructField(jsonstruct,  'diffusionModelType', 'full');
            diffusionModelType = getJsonStructField(jsonstruct, 'diffusionModelType');

            switch diffusionModelType
                
              case 'simple'
                
                jsonstruct = equalizeJsonStructField(jsonstruct, {itf, 'volumetricSurfaceArea'}, {sd, 'volumetricSurfaceArea'});
                
              case {'full', 'swelling'}
                
                jsonstruct = equalizeJsonStructField(jsonstruct, {itf, 'volumetricSurfaceArea'}, {sd, 'volumetricSurfaceArea'});

                diffusionCoefficient = getJsonStructField(jsonstruct, {sd, 'diffusionCoefficient'});
                
                if isAssigned(diffusionCoefficient)
                    % we impose that cmax in the solid diffusion model and the interface are consistent
                    paramnames = {'saturationConcentration', ...
                                  'guestStoichiometry100', ...
                                  'guestStoichiometry0'};
                    for iparam = 1 : numel(paramnames)
                        paramname = paramnames{iparam};
                        jsonstruct = equalizeJsonStructField(jsonstruct, {sd, paramname}, {itf, paramname});
                    end
                    
                end

                if strcmp(diffusionModelType, 'swelling')
                    jsonstruct = equalizeJsonStructField(jsonstruct, {sd, 'saturationConcentration'}, {itf, 'saturationConcentration'});
                end

              otherwise
                
                error('Unknown diffusionModelType %s', diffusionModelType);
                
            end

            isRootSimulationModel = getJsonStructField(jsonstruct, 'isRootSimulationModel');

            if isAssigned(isRootSimulationModel) && isRootSimulationModel 
                % only one particle in the stand-alone model
                jsonstruct = setJsonStructField(jsonstruct, {sd, 'np'}, 1);
                % For the standalone model, we set the volume fraction to one (no other component is present)
                jsonstruct = setJsonStructField(jsonstruct, {sd, 'volumeFraction'}, 1);
            end

            if jsonstruct.useLithiumPlating
                jsonstruct = equalizeJsonStructField(jsonstruct, {lp, 'volumetricSurfaceArea'}, {itf, 'volumetricSurfaceArea'});
                jsonstruct = equalizeJsonStructField(jsonstruct, {lp, 'volumeFraction'}, {sd, 'volumeFraction'});
                jsonstruct = equalizeJsonStructField(jsonstruct, {lp, 'particleRadius'}, {sd, 'particleRadius'});
            end
            
            inputparams = inputparams@ComponentInputParams(jsonstruct);
            inputparams = inputparams.setupInterface(jsonstruct);
            inputparams = inputparams.setupSolidDiffusion(jsonstruct);

            if jsonstruct.useLithiumPlating
                inputparams.LithiumPlating = LithiumPlatingLatzInputParams(pickField(jsonstruct, lp));
            end

        end

        
        function inputparams = setupInterface(inputparams, jsonstruct)

            itf = 'Interface';
            
            switch jsonstruct.SEImodel
              case {'none', 'Safari'}
                inputparams.(itf) = InterfaceInputParams(pickField(jsonstruct, itf));
              case 'Bolay'
                inputparams.(itf) = BolayInterfaceInputParams(pickField(jsonstruct, itf));
              otherwise
                error('SEImodel not recognized')
            end
            
        end
        
        function inputparams = setupSolidDiffusion(inputparams, jsonstruct)


            sd = 'SolidDiffusion';
            
            diffusionModelType = jsonstruct.diffusionModelType;
            
            switch diffusionModelType
                
              case 'simple'
                
                inputparams.(sd) = SimplifiedSolidDiffusionModelInputParams(pickField(jsonstruct, sd));
                
              case 'full'

                inputparams.(sd) = FullSolidDiffusionModelInputParams(pickField(jsonstruct, sd));

              case 'swelling'
                
                inputparams.(sd) = FullSolidDiffusionSwellingModelInputParams(pickField(jsonstruct, sd));

              otherwise
                
                error('Unknown diffusionModelType %s', diffusionModelType);
                
            end
        end
        
    end
    
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
