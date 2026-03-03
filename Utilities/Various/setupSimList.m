function [bp, simlist] = setupSimList(directory, varargin)

    opt = struct('addDirectoryName', false, ...
                 'directoryField', 'directory');
    opt = merge_options(opt, varargin{:});
    
    mrstModule add mpfa

    if exist(directory, 'dir')
        dataDirectory = directory;
    else
        dataDirectory =  fullfile(battmoDir,'Externals','mrst','core','output', directory);
    end

    assert(exist(dataDirectory, 'dir'), sprintf('Directory %s not found', dataDirectory));
    
    dataFolderObjs = dir(dataDirectory);
    dataFolderObjs = dataFolderObjs(3 : end);
    dataFolders = {};
    for ind = 1 : numel(dataFolderObjs)
        dataFolderObj = dataFolderObjs(ind);
        if dataFolderObj.isdir
            dataFolders{end + 1} = dataFolderObj.name;
        end
    end

    %% setup simlist using OutputProcessor

    bp = BatchProcessor;
    simlist = {};

    for ind = 1 : numel(dataFolders)
        dataFolder = dataFolders{ind};
        filename = fullfile(dataDirectory, dataFolder, 'input.mat');
        try
            data = load(filename);
            input = data.input;
            if opt.addDirectoryName
                input.(opt.directoryField) = dataFolder;
            end
            [bp, simlist] = bp.addElement(simlist, input);
        catch
            fprintf('no input file found for %s\n', filename);
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
