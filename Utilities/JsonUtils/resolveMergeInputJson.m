function jsonstruct = resolveMergeInputJson(jsonstruct)
% Parse the json struct, looking for every key "MergeInputs" and proceed with merging for each of them.
%    
% The object with the key "mergeInput" has the properties
% - "merge-order"   : string with one of the following
%    - "ascending"  : last inputs take precedence (default)
%    - "descending" : first inputs take precedence
% - "conflict-handling" :
%     - "warn"  : send a warning if conflict (default)
%     - "quiet" : proceed with merging quietly even in case of conflict
% - "inputs" : array of inputs to be merged
    if isstruct(jsonstruct)

        if numel(jsonstruct) == 1

            mergeInputs = getJsonStructField(jsonstruct, 'MergeInputs');

            if isAssigned(mergeInputs)

                merge_order       = getJsonStructField(mergeInputs, 'merge_order', 'ascending');
                conflict_handling = getJsonStructField(mergeInputs, 'conflict_handling', 'warn');
                inputs            = getJsonStructField(mergeInputs, 'Inputs');

                for iinput = 1 : numel(inputs)
                    jsonstructs{iinput} = resolveMergeInputJson(inputs{iinput});
                end

                if strcmp(merge_order, 'ascending')
                    jsonstructs = jsonstructs(end : -1 : 1);
                end

                warn = strcmp(conflict_handling, 'warn');
                jsonstruct = mergeJsonStructs(jsonstructs, 'warn', warn);
                
            else
                
                fds = fieldnames(jsonstruct);
                for ind = 1 : numel(fds)
                    jsonstruct.(fds{ind}) = resolveMergeInputJson(jsonstruct.(fds{ind}));
                end
                
            end

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
    

