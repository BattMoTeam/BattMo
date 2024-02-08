function jsonstruct = resolveFileInputJson(jsonstruct)
% Parse the json struct, looking for key "isFile" and replace the object there with the content of "filename" that should be given together with "isFile".
% We can give several filenames using the key "filenames" and they will be merged.
% Note that (for the moment), we do not resolve file name if they are included in json arrays (only those included in objects)


    if isstruct(jsonstruct) & numel(jsonstruct) == 1

        if ismember('isFile', fieldnames(jsonstruct)) && jsonstruct.isFile

            if ismember('filename', fieldnames(jsonstruct))
                
                % When a single entry, we convert it to cell and call function again
                jsonstruct.filenames = {jsonstruct.filename};
                jsonstruct = rmfield(jsonstruct, 'filename');
                jsonstruct = resolveFileInputJson(jsonstruct);
                
            elseif ismember('filenames', fieldnames(jsonstruct))
                
                % We resolve the first file in the list
                filename = jsonstruct.filenames{1};
                fullfilename = resolveFileName(filename);
                jsonsrc = fileread(fullfilename);
                fileJsonstruct = battMojsondecode(jsonsrc);
                fds = fieldnames(fileJsonstruct);
                for ind = 1 : numel(fds)
                    fileJsonstruct.(fds{ind}) = resolveFileInputJson(fileJsonstruct.(fds{ind}));
                end
                jsonstruct = mergeJsonStructs({fileJsonstruct, jsonstruct});

                if numel(jsonstruct.filenames) == 1
                    % If there is not filename left, we are done for this entry
                    % remove the "isFile" and "filenames" fields from the jsonstruct
                    jsonstruct = rmfield(jsonstruct, 'isFile');
                    jsonstruct = rmfield(jsonstruct, 'filenames');
                else
                    % otherwise we call the function again with the entries that are left
                    jsonstruct.filenames = jsonstruct.filenames(2 : end);
                    jsonstruct = resolveFileInputJson(jsonstruct);
                end
                
            else
                error('isFile flag is given and set to true but we are missing the filename or filenames properties so that we cannot recover the file(s)')
            end

        end

        fds = fieldnames(jsonstruct);
        for ind = 1 : numel(fds)
            jsonstruct.(fds{ind}) = resolveFileInputJson(jsonstruct.(fds{ind}));
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
