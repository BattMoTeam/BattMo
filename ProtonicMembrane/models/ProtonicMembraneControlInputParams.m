classdef ProtonicMembraneControlInputParams < InputParams
    
    properties

        controlType
        useCurrentDensity % If true, then we use the value of currentDensity. Otherwise, a value for I the total current should
                          % be given.
        currentDensity
        I
        
        area % if current density is used, we need the area. This value will typically be assigned by the grid constructor
        
    end
    
    methods
        
        function inputparams = ProtonicMembraneControlInputParams(jsonstruct)

            jsonstruct  = setDefaultStructField(jsonstruct, {'useCurrentDensity'}, false);
            jsonstruct  = setDefaultStructField(jsonstruct, {'I'}, 0);
            inputparams = inputparams@InputParams(jsonstruct);
            
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
