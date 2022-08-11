function res = decomposeVarName(varname, givenvarname)
% check if varname matches with givenvarname.
%
% varname and given varname matches but with different indices, return structure res
%    
%   res.found                : if there is match (same name but different indices)
%   res.intersecVarname      : VarName with the common indices
%   res.reminderGivenVarname : VarName with the indices that are in givenvarname and not in varname
%   res.reminderVarname      : VarName with the indices that are in varname but not in givenvarname
%

    if hasSameName(varname, givenvarname)
        
        ind = varname.index;
        gind = givenvarname.index;
        
        if iscolumn(ind) & iscolumn(gind)
            res = struct('found'               , true   , ...
                        'intersecVarname'     , varname, ...
                        'reminderVarname'     , []     , ...
                        'reminderGivenVarname', []);
            return
        end

        if iscolumn(ind)
            ind = (1 : varname.dim)';
        end
        
        if iscolumn(gind)
            gind = (1 : givenvarname.dim)';
        end
        
        [lia, locb] = ismember(ind, gind);

        if all(lia == 0)

            res = struct('found'               , false   , ...
                        'intersecVarname'     , [], ...
                        'reminderVarname'     , varname     , ...
                        'reminderGivenVarname', givenvarname);
            return

        else

            % Initialize the structure
            intersecVarname      = varname;
            reminderVarname      = varname;
            reminderGivenVarname = varname;
            
        end
        
        intersecVarname.index = ind(lia);
        if all(lia)
            reminderVarname = [];
        else
            reminderVarname.index = ind(~lia);
        end
        
        [lia, locb] = ismember(gind, ind);        
        if all(lia)
            reminderGivenVarname = [];
        else
            reminderGivenVarname.index = gind(~lia);
        end
        
        res = struct('found'               , true           , ...
                    'intersecVarname'     , intersecVarname, ...
                    'reminderVarname'     , reminderVarname, ...
                    'reminderGivenVarname', reminderGivenVarname);
        
    else

        res = struct('found'               , false  , ...
                    'intersecVarname'     , []     , ...
                    'reminderVarname'     , varname, ...
                    'reminderGivenVarname', givenvarname);
        
    end        



end

function is = iscolumn(val)
    
    if ischar(val) && strcmp(val, ':')
        is = true;
    else
        is = false;
    end
    
end

function isequal = hasSameName(varname1, varname2)
% check only on namespace and name and not on index
    
    namespace1 = varname1.namespace;
    name1      = varname1.name;
    
    namespace2 = varname2.namespace;
    name2      = varname2.name;
    
    isequal = true;
    
    if ~strcmp(name1, name2)
        isequal = false;
        return
    end

    if numel(namespace1) ~= numel(namespace2)
        isequal = false;
        return
    end
    
    for i = 1 : numel(namespace1)
        if ~strcmp(namespace1{i}, namespace2{i})
            isequal = false;
            return
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
