classdef SeaWaterBatteryInputParams < ComponentInputParams

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

        function inputparams = SeaWaterBatteryInputParams(jsonstruct)

            jsonstruct = setDefaultJsonStructField(jsonstruct, 'include_precipitation', true);

            inputparams = inputparams@ComponentInputParams(jsonstruct);
            
            inputparams.Cathode               = HydrogenElectrodeInputParams(jsonstruct.Cathode);
            inputparams.CathodeActiveMaterial = HydrogenActiveMaterialInputParams(jsonstruct.CathodeActiveMaterial);
            inputparams.Anode                 = SeaWaterElectrodeInputParams(jsonstruct.Anode);
            inputparams.AnodeActiveMaterial   = MagnesiumActiveMaterialInputParams(jsonstruct.AnodeActiveMaterial);
            
            if inputparams.include_precipitation
                inputparams.Electrolyte = SeaWaterElectrolyteInputParams(jsonstruct.Electrolyte);
            else
                inputparams.Electrolyte = SeaWaterElectrolyteNoPrecipitationInputParams(jsonstruct.Electrolyte);
            end
            
            inputparams.couplingTerms = {};

        end

        function inputparams = validateInputParams(inputparams)

            if inputparams.include_precipitation
                assert(isa(inputparams.Electrolyte, 'SeaWaterElectrolyteInputParams'), 'The input class SeaWaterElectrolyteNoPrecipitationInputParams is used. The input parameters for the electrolyte may include precipitation data.');
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
