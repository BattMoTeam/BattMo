classdef CoatingInputParams < ElectronicComponentInputParams
%
% Input parameter class for :code:`Coating` model
%
    properties

        %% Sub-Models

        ActiveMaterial
        Binder
        ConductingAdditive

        % The two following models are instantiated only when active_material_type == 'composite' and, in this case,
        % ActiveMaterial model will remain empty. If active_material_type == 'default', then the two models remains empty
        ActiveMaterial1
        ActiveMaterial2

        %% Standard parameters

        effectiveDensity     % the mass density of the material (symbol: rho). Important : the density is computed in wet or calendared state (meaning including the volume of the pores)
        bruggemanCoefficient % the Bruggeman coefficient for effective transport in porous media (symbol: beta)
        active_material_type % 'default' (only one particle type) or 'composite' (two different particles)

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

            if isempty(inputparams.active_material_type)
                inputparams.active_material_type = 'default';
            end

            switch inputparams.active_material_type
              case 'default'

                am = 'ActiveMaterial';
                inputparams.(am) = ActiveMaterialInputParams(jsonstruct.(am));

              case 'composite'

                am1 = 'ActiveMaterial1';
                am2 = 'ActiveMaterial2';
                inputparams.(am1) = ActiveMaterialInputParams(jsonstruct.(am1));
                inputparams.(am2) = ActiveMaterialInputParams(jsonstruct.(am2));

              otherwise
                error('active_material_type not recognized');
            end
            inputparams.Binder             = BinderInputParams(pick('Binder'));
            inputparams.ConductingAdditive = ConductingAdditiveInputParams(pick('ConductingAdditive'));

            inputparams = inputparams.validateInputParams();

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
