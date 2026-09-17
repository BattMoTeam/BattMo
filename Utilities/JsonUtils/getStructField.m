function value = getStructField(jsonstruct, fieldnamelist, defaultValue)

    if ischar(fieldnamelist)
        % handle case where fieldnamelist is just a char
        fieldnamelist = {fieldnamelist};
        if nargin > 2
            value = getStructField(jsonstruct, fieldnamelist, defaultValue);
        else
            value = getStructField(jsonstruct, fieldnamelist);
        end
        
        return
    end

    fieldname = fieldnamelist{1};

    if numel(fieldnamelist) > 1

        fieldnamelist = fieldnamelist(2 : end);
        getValue = false;
        
    else
        
        getValue = true;
        
    end

    if isempty(jsonstruct) || (~isfield(jsonstruct, fieldname) && ~isprop(jsonstruct, fieldname))

        value = UnAssigned();

    else

        if getValue

            value = jsonstruct.(fieldname);

        else

            value = getStructField(jsonstruct.(fieldname), fieldnamelist);

        end
        
    end

    if isUnAssigned(value) && nargin > 2
        
        value = defaultValue;
        
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
