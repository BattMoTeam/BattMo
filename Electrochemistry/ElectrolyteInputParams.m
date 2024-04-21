classdef ElectrolyteInputParams < ComponentInputParams
%
% Input parameter class for :code:`Electrolyte` model
%
    properties

        species % Structure given properties of the main ion (Lithium)
        
        %% Standard parameters

        density              % the mass density of the material (symbol: rho)
        ionicConductivity    % a function to determine the ionic conductivity of the electrolyte under given conditions (symbol: kappa)
        diffusionCoefficient % a function to determine the diffusion coefficient of a molecule in the electrolyte under given conditions (symbol: D)        
        bruggemanCoefficient % the coefficient for determining effective transport parameters in porous media (symbol: beta)

        thermalConductivity  % Intrinsic Thermal conductivity of the electrolyte
        specificHeatCapacity % Specific Heat capacity of the electrolyte

        nominalEthyleneCarbonateConcentration % only used if a SEI layer is included in the model
        
        %% Advanced parameters

        volumeFraction
        effectiveThermalConductivity    % (account for volume fraction)
        effectiveVolumetricHeatCapacity % (account for volume fraction and density)
        
        use_thermal
        
    end

    methods

        function inputparams = ElectrolyteInputParams(jsonstruct)

            inputparams = inputparams@ComponentInputParams(jsonstruct);

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
