classdef TimeControlModelInputParams < ControlModelInputParams

    properties

        value % function setup (see Utilities/JsonSchemas/Function.schema.json) that returns the value for the control,
              % either current or voltage, depending on the returned value of the type function
        type % function setup (see Utilities/JsonSchemas/Function.schema.json) that returns the value for the control
             % type,
             % - 1 for current
             % - 2 for voltage

        % Advanced parameters

        tolerance = 1e-4 % tolerance to skip timesteps (in second)

    end

    methods
        
        function inputparams = TimeControlModelInputParams(jsonstruct);

            inputparams =  inputparams@ControlModelInputParams(jsonstruct);
            
        end
        
    end
    
end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
