classdef GenericControlModelInputParams < ControlModelInputParams
%
    properties

        controlsteps
        capacity % value is used if CRate or DRate are given

        addDescription % boolean, if true, a description of the step is added to each state. default is false
        description % description string
        
    end
    
    
    methods

        function inputparams = GenericControlModelInputParams(jsonstruct)

            jsonstruct = setDefaultJsonStructField(jsonstruct, {'addDescription'}, false);
            jsonstruct = setDefaultJsonStructField(jsonstruct, {'description'}, '');
            
            % addon to compensate for matlab json parsing choice
            if isstruct(jsonstruct.controlsteps)
                jsonstruct.controlsteps = {jsonstruct.controlsteps};
            end
            inputparams = inputparams@ControlModelInputParams(jsonstruct);
            
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
