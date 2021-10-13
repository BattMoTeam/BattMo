function jsonstruct = resolveFileInputJson(jsonstruct)
    
    fileroot = batmoDir();
    
    if isstruct(jsonstruct)
        fds = fieldnames(jsonstruct);
        if ismember('isFile', fds)
            filename = jsonstruct.filename;
            fullfilename = fullfile(fileroot, filename);
            jsonsrc = fileread(fullfilename);
            parsedjson = jsondecode(jsonsrc);
            jsonstruct = parsedjson;
        end
        fds = fieldnames(jsonstruct);
        for ind = 1 : numel(fds)
            jsonstruct.(fds{ind}) = resolveFileInputJson(jsonstruct.(fds{ind}));
        end
    end
    
end
