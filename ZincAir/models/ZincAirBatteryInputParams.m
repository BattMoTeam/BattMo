classdef ZincAirBatteryInputParams < ComponentInputParams

    properties

        Cathode
        CathodeActiveMaterial
        Anode
        AnodeActiveMaterial
        Electrolyte
        
        couplingTerms    % list of coupling terms (each element is instance of couplingTerm class)

        T % Temperature is given as input (for the moment)

        include_precipitation
    end
    
    methods

        function inputparams = ZincAirBatteryInputParams(jsonstruct)

            inputparams = inputparams@ComponentInputParams(jsonstruct);
            
            inputparams.Cathode               = OxygenElectrodeInputParams(jsonstruct.Cathode);
            inputparams.CathodeActiveMaterial = ZincAirActiveMaterialInputParams(jsonstruct.CathodeActiveMaterial);
            inputparams.Anode                 = ZincAirElectrodeInputParams(jsonstruct.Anode);
            inputparams.AnodeActiveMaterial   = ZincActiveMaterialInputParams(jsonstruct.AnodeActiveMaterial);
            if inputparams.include_precipitation
                inputparams.Electrolyte = ZincAirElectrolyteInputParams(jsonstruct.Electrolyte);
            else
                inputparams.Electrolyte = ZincAirElectrolyteInputParams(jsonstruct.Electrolyte);
            end
            
            inputparams.couplingTerms = {};

        end

        function inputparams = validateInputParams(inputparams)

            if inputparams.include_precipitation
                assert(isa(inputparams.Electrolyte, 'ZincAirElectrolyteInputParams'), 'The input class ZincAirElectrolyteInputParams is used. The input parameters for the electrolyte may include precipitation data.');
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
