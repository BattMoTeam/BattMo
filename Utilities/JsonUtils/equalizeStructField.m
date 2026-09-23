function [jsonstruct, bothUnAssigned] = equalizeStructField(jsonstruct, fieldnamelist1, fieldnamelist2, varargin)

    value1 = getStructField(jsonstruct, fieldnamelist1);
    value2 = getStructField(jsonstruct, fieldnamelist2);

    bothUnAssigned = false;
    
    if isUnAssigned(value1)

        if isUnAssigned(value2)
            
            bothUnAssigned = true;
            return
            
        else
            
            jsonstruct = setStructField(jsonstruct, fieldnamelist1, value2);

        end

    else

        if isUnAssigned(value2)
            
            jsonstruct = setStructField(jsonstruct, fieldnamelist2, value1);

        else

            if isequal(value1, value2)

                return

            else

                opt = struct('force', 'false', ...
                             'warn', true);
                opt = merge_options(opt, varargin{:});

                if opt.force
                    errorMessage = sprintf('Different values are given for the fields jsonstruct.%s and jsonstruct.%s. We do not know which one to choose...', ...
                                           getPrintableName(fieldnamelist1)                                                                                  , ...
                                           getPrintableName(fieldnamelist2));
                    jsonstruct = setStructField(jsonstruct, fieldnamelist2, value1, 'handleMisMatch', 'error', 'errorMessage', errorMessage);
                    if opt.warn
                        fprintf('Fist value given in equalizeStructField is taken\n');
                    end
                else
                    error('mismatch in equalizeStructField');
                end
                
            end
        end

    end
      
end


function namestr = getPrintableName(fieldnamelist)

    if ischar(fieldnamelist)

        namestr = getPrintableName({fieldnamelist})
        return
        
    end

    namestr = strjoin(fieldnamelist, '.')

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
