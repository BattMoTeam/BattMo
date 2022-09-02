classdef ControlModelInputParams < InputParams

    
    properties

        controlPolicy
        CRate
        
    end
    
    methods

        function paramobj = ControlModelInputParams(jsonstruct);
            
            paramobj = paramobj@InputParams(jsonstruct);
            
        end
        
        function paramobj = set.controlPolicy(paramobj, controlPolicy)
            switch controlPolicy
              case 'IEswitch'
                % ok in any case
              case 'CCCV'
                assert(isa(paramobj, 'CcCvControlModelInputParams'), 'The model is not a CcCvControlModelInputParams class')
              case 'CV'
                assert(isa(paramobj, 'CvControlModelInputParams'), 'The model is not a CvControlModelInputParams class')
              otherwise
                error('controlPolicy not recognized');
            end
            paramobj.controlPolicy = controlPolicy;                
        end
        
    end
    
end



%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
