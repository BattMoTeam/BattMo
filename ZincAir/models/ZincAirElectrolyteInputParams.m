classdef ZincAirElectrolyteInputParams < ElectronicComponentInputParams

    properties

        K % Reaction rates for the following reactions (depends on system)
        
        kappa % conductivity
        
        species % cell array for species 
                % Each cell is a struct with fields
                % - name (string)
                % - z (only for solutes)
                % - lambda0 (only for solutes)
                % - D : diffusion coefficients

        solutes
        solids
        logspecies
        
        quasiparticles % cell array of struct with fields
                       % - name (string)
                       % - composition : array of struct with field
                       %       - name : name of species
                       %       - coef : coefficient in the quasiParticle decomposition

        
    end
    
    methods
        
        
        function inputparams = ZincAirElectrolyteInputParams(jsonstruct)
            
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
