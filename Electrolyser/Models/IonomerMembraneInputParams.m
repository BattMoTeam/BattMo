classdef IonomerMembraneInputParams < ElectronicComponentInputParams
    
    properties
        
        volumeFraction
        
        H2O % with fields
        % H2O.c0 : Reference concentration
        % H2O.D : diffusion coefficient for water
        % H2O.MW : Molar mass (needed for function groupHydration which is only needed in setup of initial condition and not for assembly)

        
        OH % with fields
        % OH.xi : OH occupation
        % OH.z : Charge number
        % OH.t : Transference number

        cT % Total concentration of charged group (one scalar value)
        
        V % molar volume (needed for function groupHydration which is only needed in setup of initial condition and not for assembly)

        tortuosity % cell-valued coefficient
        
    end
    
    methods
        
        function paramobj = IonomerMembraneInputParams(jsonstruct)
            
            paramobj = paramobj@ElectronicComponentInputParams(jsonstruct);
            
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
