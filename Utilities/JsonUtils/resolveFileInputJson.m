function jsonstruct = resolveFileInputJson(jsonstruct)
% Parse the json struct, looking for key "isFile" and replace the object there with the content of "filename" that
% should be given together with "isFile".


    if isstruct(jsonstruct)

        if numel(jsonstruct) == 1

            isFile = getJsonStructField(jsonstruct, 'isFile');

            if isAssigned(isFile)
                
                filename  = getJsonStructField(jsonstruct, 'filename');
                filenames = getJsonStructField(jsonstruct, 'filenames');
                
                if isAssigned(filename)
                    
                    fullfilename = resolveFileName(filename);
                    jsonsrc = fileread(fullfilename);
                    fileJsonstruct = battMojsondecode(jsonsrc);
                    fds = fieldnames(fileJsonstruct);
                    for ind = 1 : numel(fds)
                        fileJsonstruct.(fds{ind}) = resolveFileInputJson(fileJsonstruct.(fds{ind}));
                    end
                    jsonstruct = fileJsonstruct;

                elseif isAssigned(filenames)
                    
                    filenames = jsonstruct.filenames;

                    jsonstruct = removeJsonStructField(jsonstruct, {'isFile'}, 'handleMissing', 'quiet');
                    jsonstruct = removeJsonStructField(jsonstruct, {'filenames'}, 'handleMissing', 'quiet');
                    MergeInputs = struct('merge_order', 'ascending', ...
                                         'conflict_handling', 'warn');
                    Inputs = {};
                    for ifile = 1 : numel(filenames)
                        clear fileinput
                        fileinput.isFile = true;
                        fileinput.filename = filenames{ifile};
                        Inputs{ifile} = fileinput;
                    end
                    MergeInputs.Inputs = Inputs;
                    jsonstruct.MergeInputs = MergeInputs;
                    jsonstruct = resolveFileInputJson(jsonstruct);
                    return

                else
                    error('missing filename or filenames keyword');
                end
                

            end

            fds = fieldnames(jsonstruct);
            for ind = 1 : numel(fds)
                jsonstruct.(fds{ind}) = resolveFileInputJson(jsonstruct.(fds{ind}));
            end
            
        else
            
            % we have an array. We call resolveFileInputJson on each element
            for i_jsonstruct = 1 : numel(jsonstruct)
                jsonstruct(i_jsonstruct) = resolveFileInputJson(jsonstruct(i_jsonstruct));
            end
           
        end

    elseif iscell(jsonstruct)

        % we have cells. We call resolveFileInputJson on each element
        for i_jsonstruct = 1 : numel(jsonstruct)
            jsonstruct{i_jsonstruct} = resolveFileInputJson(jsonstruct{i_jsonstruct});
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
