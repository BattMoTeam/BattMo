classdef InputParams 
    
    methods

        function paramobj = InputParams(jsonstruct)
            paramobj = assignJsonParams(paramobj, jsonstruct);
        end
        
        function paramobj = setParam(paramobj, names, val)
        % Note : same syntax as setProp in BaseModel
            if iscell(names) & (numel(names) > 1)
                name = names{1};
                names = names(2 : end);
                paramobj.(name) = setParam(paramobj.(name), names, val);
            elseif iscell(names) & (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    paramobj{name} = val;
                else
                    paramobj.(name) = val;
                end
            else
                error('format not recognized');
            end
        end

        function var = getParam(paramobj, names)
        % Note : same syntax as getProp in BaseModel
            if iscell(names) && (numel(names) > 1)
                name = names{1};
                names = names(2 : end);
                var = getParam(paramobj.(name), names);
            elseif iscell(names) & (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    var = paramobj{name};
                else
                    var = paramobj.(name);
                end
            else
                error('format not recognized');
            end
        end

        function paramobj = validateInputParams(paramobj)

        % Default automatic behaviour is that all the properties of paramobj that belong to class InputParams get
        % validated.
            
            paramobjFds = propertynames(paramobj);

            for ind = 1 : numel(paramobjFds)
        
                fd = paramobjFds{ind};
        
                if isa(paramobj.(fd), 'InputParams')
                    paramobj.(fd) = paramobj.(fd).validateInputParams();
                end
                
            end
    
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
