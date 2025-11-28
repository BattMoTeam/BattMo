classdef InterfaceInputParams < InputParams 
%
% Input parameter class for the interface model
%    
    properties

        G  % Grid

        saturationConcentration      % the saturation concentration of the guest molecule in the host material
        numberOfElectronsTransferred % stoichiometric number of electrons transferred in the electrochemical reaction
        volumetricSurfaceArea        % surface area of the active material - electrolyte interface per volume of electrode
        activationEnergyOfReaction   % the activation energy of the electrochemical reaction
        reactionRateConstant         % the reaction rate constant of the electrochemical reaction

        % The exchange current at the electrode-electrolyte interface under equilibrium conditions.
        % Tf empty, the default expression using the reaction rate constant is used, see method
        % Interface.updateReactionRateCoefficient. The function is given as a struct with the fields:
        %   - type = {"function", "constant"} % if "constant" is selected, we use the reactionRateConstant value
        %   - functionName :  matlab function name (should be available in path)
        %   - argumentList = ["cElectrodeSurface", "cmax"]
        exchangeCurrentDensity
        
        guestStoichiometry100 % the ratio of the concentration of the guest molecule to the saturation concentration
                              % of the guest molecule in a phase at a cell voltage that is defined as 100% SOC
        guestStoichiometry0   % the ratio of the concentration of the guest molecule to the saturation concentration
                              % of the guest molecule in a phase at a cell voltage that is defined as 0% SOC
        density               % the mass density of the active material

        % A function to determine the open-circuit potential of the electrode under given conditions
        % See schema `Utilities/JsonSchemas/Function.schema.json` for a complete description of the function interface
        openCircuitPotential

        includeEntropyChange     % flag which determines if entropy change should be computed and included

        % A function to determine the entropy change
        % See schema `Utilities/JsonSchemas/Function.schema.json` for a complete description of the function interface        
        entropyChange

        referenceTemperature % Used to compute effective OCP from reference OCP and entropy change
        
        chargeTransferCoefficient % the charge transfer coefficient that enters in the Butler-Volmer equation (symbol: alpha)

        %% SEI model
        SEImodel % string defining interface type. Can take value
                  % - 'none' (default)
                  % - 'Safari'
                  % - 'Bolay'
        
        %% Double layer capacity
        useDoubleLayerCapacity % if true, add double layer capacity (default is false)
        doubleLayerCapacitance % Value of electric double layer capacitance / Fm^-2      
    end
    
    methods
        
        function inputparams = InterfaceInputParams(jsonstruct)

            if isUnAssigned(jsonstruct, 'entropyChange')
                jsonstruct = setJsonStructField(jsonstruct, 'includeEntropyChange', false);
            else
                jsonstruct = setDefaultJsonStructField(jsonstruct, 'includeEntropyChange', true);
            end
            
            inputparams = inputparams@InputParams(jsonstruct);
            
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
