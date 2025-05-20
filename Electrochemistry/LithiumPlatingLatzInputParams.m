classdef LithiumPlatingLatzInputParams < InputParams
%
% Input parameter class 
% 
    properties

        alphaPl   
        alphaStr  
        alphaChInt

        kPl   
        kChInt

        nPl0
        nPlLimit

        SEIFraction
        MSEI       
        rhoSEI     
        deltaSEI0  
        sigmaSEI   

        useSEI
        
    end

    methods

        function inputparams = LithiumPlatingLatzInputParams(jsonstruct)

            jsonstruct = setDefaultJsonStructField(jsonstruct, 'alphaPl'    , 0.3);
            jsonstruct = setDefaultJsonStructField(jsonstruct, 'alphaStr'   , 0.7);
            jsonstruct = setDefaultJsonStructField(jsonstruct, 'alphaChInt' , 0.5);
            jsonstruct = setDefaultJsonStructField(jsonstruct, 'kPl'        , 1e-10);
            jsonstruct = setDefaultJsonStructField(jsonstruct, 'kChInt'     , 1e-12);
            jsonstruct = setDefaultJsonStructField(jsonstruct, 'nPl0'       , 1.173e-23);
            jsonstruct = setDefaultJsonStructField(jsonstruct, 'nPlLimit'   , 1.173e-17);
            jsonstruct = setDefaultJsonStructField(jsonstruct, 'SEIFraction', 0.05);
            jsonstruct = setDefaultJsonStructField(jsonstruct, 'MSEI'       , 0.162);
            jsonstruct = setDefaultJsonStructField(jsonstruct, 'rhoSEI'     , 1690);
            jsonstruct = setDefaultJsonStructField(jsonstruct, 'deltaSEI0'  , 1e-9);
            jsonstruct = setDefaultJsonStructField(jsonstruct, 'sigmaSEI'   , 5e-6);
            jsonstruct = setDefaultJsonStructField(jsonstruct, 'useSEI'     , false);
            
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
