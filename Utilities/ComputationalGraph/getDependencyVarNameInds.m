function [depvarnameinds, propfuncinds, propvarnameinds, staticinds] = getDependencyVarNameInds(varnameinds, A)
    
    [propfuncinds, propvarnameinds, staticinds] = getDependencyVarNameIndsRecursive(varnameinds, [], [], [], A);

    depvarnameinds = [staticinds; propvarnameinds];
    
end

function [propfuncinds, propvarnameinds, staticinds] = getDependencyVarNameIndsRecursive(varnameinds, propfuncinds, propvarnameinds, staticinds, A)
% for given list of index varnameinds, returns a list with index of nodes

    depvarnameinds = [];

    for ivar = 1 : numel(varnameinds)

        varnameind = varnameinds(ivar);
        c = A(:, varnameind);

        if all(c == 0)
            %  no dependency
            staticinds = vertcat(staticinds, varnameind);
        else
            % the non-zero entry in the column is is unique and give the index of the property function
            
            propfuncind = max(c);

            depvarnameinds  = vertcat(depvarnameinds, find(c));
            propfuncinds    = vertcat(propfuncind, propfuncinds);
            propvarnameinds = vertcat(varnameind, propvarnameinds);
            
        end
        
    end        
    
    if ~isempty(depvarnameinds)

        varnameinds = depvarnameinds;
        [propfuncinds, propvarnameinds, staticinds] = getDependencyVarNameIndsRecursive(varnameinds, propfuncinds, propvarnameinds, staticinds, A);
        
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
