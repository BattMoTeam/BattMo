function jsonstruct = setStructField(jsonstruct, fieldnamelist, value, varargin)
    % in varargin, the key 'handleMisMatch' can take following values
    % - 'error'   : returns error if value already set and does not match with the new given one (default)
    % - 'quiet'   : does not warn about case above
    % - 'warning' : warn but no error

    if ischar(fieldnamelist)
        % handle case wher fieldnamelist is just a char
        fieldnamelist = {fieldnamelist};
        jsonstruct = setStructField(jsonstruct, fieldnamelist, value, varargin{:});
        return
    end

    fieldname = fieldnamelist{1};

    if numel(fieldnamelist) > 1
        fieldnamelist = fieldnamelist(2 : end);
        setValue = false;
    else
        setValue = true;
    end

    if isAssigned(jsonstruct, fieldname)

        if setValue

            currentValue = jsonstruct.(fieldname);

            equalValue = isequal(currentValue, value);
            
            if equalValue
                
                % do nothing. Value was already set to given value
                return
                
            else
                
                opt = struct('handleMisMatch', 'error', ...
                             'errorMessage', []);
                opt = merge_options(opt, varargin{:});
                
                switch opt.handleMisMatch

                  case 'quiet'

                    jsonstruct.(fieldname) = value;
                    
                  case 'warning'

                    fprintf('mismatch values in assignment of %s. We use the given value\n', fieldname)
                    jsonstruct.(fieldname) = value;
                    
                  case 'error'

                    if isempty(opt.errorMessage)
                        errorMessage = sprintf('mismatch values in assignment of %s. We use the given value\n', fieldname);
                    else
                        errorMessage = opt.errorMessage;
                    end
                    error(errorMessage);

                  otherwise

                    error('handleMisMatch not recognized');
                    
                end

            end
        else

            jsonstruct.(fieldname) = setStructField(jsonstruct.(fieldname), fieldnamelist, value, varargin{:});

        end

    else

        if setValue

            jsonstruct.(fieldname) = value;

        else

            jsonstruct.(fieldname) = setStructField([], fieldnamelist, value, varargin{:});
            
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
