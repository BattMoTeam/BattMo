function jsonstruct = resolveFileInputJson(jsonstruct)
% Note that (for the moment), we do not resolve file name if they are included in json arrays (only those included in objects)
    
    fileroot = battmoDir();
    
    if isstruct(jsonstruct) & numel(jsonstruct) == 1

        if ismember('isFile', fieldnames(jsonstruct))
            %% TODO : check if it isFile is true
            filename = jsonstruct.filename;
            fullfilename = fullfile(fileroot, filename);
            jsonsrc = fileread(fullfilename);
            fileJsonstruct = jsondecode(jsonsrc);
            
            % remove the "isFile" and "filename" fields from the jsonstruct
            jsonstruct = rmfield(jsonstruct, 'isFile');
            jsonstruct = rmfield(jsonstruct, 'filename');
            
            fds = fieldnames(fileJsonstruct);
            for ind = 1 : numel(fds)
                fileJsonstruct.(fds{ind}) = resolveFileInputJson(fileJsonstruct.(fds{ind}));
            end

            jsonstruct = mergeJsonStructs({fileJsonstruct, jsonstruct});

            
        end
        
        fds = fieldnames(jsonstruct);
        for ind = 1 : numel(fds)
            jsonstruct.(fds{ind}) = resolveFileInputJson(jsonstruct.(fds{ind}));
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
