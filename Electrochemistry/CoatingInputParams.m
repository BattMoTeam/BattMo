classdef CoatingInputParams < ElectronicComponentInputParams
%
% Input parameter class for :code:`Coating` model
%
    properties

        %% Sub-Models

        ActiveMaterial
        Binder
        ConductingAdditive

        % The two following models are instantiated only when activeMaterialModelSetup.composite is true and, in this case,
        % ActiveMaterial model will remain empty. If activeMaterialModelSetup.composite is false, then the two models remains empty
        ActiveMaterial1
        ActiveMaterial2

        %% Standard parameters

        effectiveDensity % the mass density of the material (symbol: rho). Important : the density is computed in wet or
                         % calendared state (meaning including the volume of the pores)

        bruggemanCoefficient % the Bruggeman coefficient for effective transport in porous media (symbol: beta)
        
        activeMaterialModelSetup % instance of ActiveMaterialModelSetupInputParams. Contains fields
                              % - 'composite' : boolean (default is false)
                              % - 'SEImodel' : string with one of
                              %                 "none" (default)
                              %                 "Safari"
                              %                 "Bolay"

        %% Advanced parameters

        volumeFractions
        volumeFraction
        thermalConductivity             % (if not given computed from the subcomponents)
        specificHeatCapacity            % (if not given computed from the subcomponents)
        effectiveThermalConductivity    % (account for volume fraction)
        effectiveVolumetricHeatCapacity % (account for volume fraction and density)

        %% External coupling parameters

        externalCouplingTerm % structure to describe external coupling (used in absence of current collector)

    end

    methods

        function inputparams = CoatingInputParams(jsonstruct)

            inputparams = inputparams@ElectronicComponentInputParams(jsonstruct);

            pick = @(fd) pickField(jsonstruct, fd);

            if isfield(jsonstruct, 'activeMaterialModelSetup')
                inputparams.activeMaterialModelSetup  = ActiveMaterialModelSetupInputParams(jsonstruct.activeMaterialModelSetup);
            else
                inputparams.activeMaterialModelSetup  = ActiveMaterialModelSetupInputParams([]);
            end
            
            if inputparams.activeMaterialModelSetup.composite
                
                am1 = 'ActiveMaterial1';
                am2 = 'ActiveMaterial2';
                inputparams.(am1) = ActiveMaterialInputParams(jsonstruct.(am1));
                inputparams.(am2) = ActiveMaterialInputParams(jsonstruct.(am2));

            else
                
                am = 'ActiveMaterial';

                switch inputparams.activeMaterialModelSetup.SEImodel

                  case {'none', 'Balay'}

                    inputparams.(am) = ActiveMaterialInputParams(jsonstruct.(am));
                    
                  case 'Safari'

                    inputparams.(am) = SEIActiveMaterialInputParams(jsonstruct.(am));

                  otherwise
                    
                    error('active material modelSEI layer model not recognized');
                    
                end
            end
            
            inputparams.Binder             = BinderInputParams(pick('Binder'));
            inputparams.ConductingAdditive = ConductingAdditiveInputParams(pick('ConductingAdditive'));

        end

        function inputparams = validateInputParams(inputparams)

            inputparams = validateInputParams@InputParams(inputparams);
            inputparams = mergeParameters(inputparams, {{'ActiveMaterial', 'SEImodel'}, {'ActiveMaterial', 'Interface', 'SEImodel'}});
            

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
