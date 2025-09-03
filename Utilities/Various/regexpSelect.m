function inds = regexpSelect(names, name)
% implementation of regexp where whitespace matches every thing and name can be a cell in which case the match of each element will be included

    if iscell(name)
        inds = [];
        for iname = 1 : numel(name)
            addedinds = regexpSelect(names, name{iname});
            inds = vertcat(inds, addedinds);
        end
        inds = unique(inds);
    else
        name = replace(name, '{', '\{');
        name = replace(name, '}', '\}');
        name = regexprep(name, ' +', '.*');
        inds = regexp(names, name, 'once');
        inds = cellfun(@(x) ~isempty(x), inds);
        inds = find(inds);
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
