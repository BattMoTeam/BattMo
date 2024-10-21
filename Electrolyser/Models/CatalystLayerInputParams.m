classdef CatalystLayerInputParams < ComponentInputParams
    
    properties

        referenceExchangeCurrentDensity    % Exchange current density
        standardEquilibriumPotential    % Standard equilibrium potential
        referencePotential  % Reference potential

        species % species struct with field
        % - OH.chargeNumber  : Charge number
        % - OH.referenceConcentration : OH reference concentration

        numberOfElectronsTransferred % Number of electron transfer
        
        chargeTransferCoefficient                  % coefficient in the exponent in Butler-Volmer equation [-]
        ionomerFractionArea                  % Fraction of specific area that is coversed with ionomer [-]
        referenceVolumetricSurfaceArea % Volumetric surface area [m^ -1]

        tortuosity % Tortuosity [-]

        include_dissolution % True if dissolution model is included (if not given, it is set to false)

        DissolutionModel
        
    end
    
    methods
        
        function inputparams = CatalystLayerInputParams(jsonstruct)

            inputparams = inputparams@ComponentInputParams(jsonstruct);

            if isempty(inputparams.include_dissolution)
                inputparams.include_dissolution = false;
            elseif inputparams.include_dissolution
                inputparams.DissolutionModel = DissolutionModelInputParams(jsonstruct.DissolutionModel);
            end

        end

        function inputparams = validateInputParams(inputparams)

            inputparams = validateInputParams@ComponentInputParams(inputparams);
            
            dm = 'DissolutionModel';

            if inputparams.include_dissolution
                inputparams = mergeParameters(inputparams, {{'referenceVolumetricSurfaceArea'}, {dm, 'referenceVolumetricSurfaceArea'}});
                inputparams.(dm) = inputparams.(dm).validateInputParams();
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
