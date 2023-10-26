function [bp, simlist] = setupSimList(directory)

    mrstModule add   mpfa

    dataDirectory =  fullfile('/home/xavier/Matlab/Projects/battmo/mrst/mrst-core/output/', directory);
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
            [bp, simlist] = bp.addElement(simlist, input);
        catch
            fprintf('no input file found for %s\n', filename);
        end
    end

end