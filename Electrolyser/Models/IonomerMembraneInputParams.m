classdef IonomerMembraneInputParams < ElectronicComponentInputParams
    
    properties
        
        volumeFraction
        
        H2O % with fields
        % H2O.referenceConcentration : Reference concentration
        % H2O.diffusionCoefficient : diffusion coefficient for water
        % H2O.molecularWeight : Molar mass (needed for function groupHydration which is only needed in setup of initial condition and not for assembly)

        
        OH % with fields
        % OH.occupationNumber : OH occupation
        % OH.chargeNumber : Charge number
        % OH.transference : Transference number

        totalChargedGroupConcentration % Total concentration of charged group (one scalar value)
        
        molarVolume % molar volume (needed for function groupHydration which is only needed in setup of initial condition and not for assembly)

        tortuosity % cell-valued coefficient
        
    end
    
    methods
        
        function inputparams = IonomerMembraneInputParams(jsonstruct)
            
            inputparams = inputparams@ElectronicComponentInputParams(jsonstruct);
            
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
