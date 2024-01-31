function data = class2data(obj)
% Convert Matlab object into struct which can be then read in Julia
    
    if(isstruct(obj) || isobject(obj))
        myfields = fields(obj);
        if(numel(obj)>1)
            assert(isstruct(obj))
            data = obj;
            for k=1:numel(obj)  
                for i=1:numel(myfields)
                    data(k).(myfields{i}) = class2data(obj(k).(myfields{i}));            
                end
            end
        else
            data = struct();
            for i=1:numel(myfields)
                data.(myfields{i}) = class2data(obj.(myfields{i}));            
            end
        end
    else
        if(isnumeric(obj) || ischar(obj) )
            data = obj;
        elseif (iscell(obj))
            vals = obj;
            newvals ={};
            for i=1:numel(vals)
                newvals{i} = class2data(vals{i});
            end
            data = newvals;
        else
            data = [];
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
