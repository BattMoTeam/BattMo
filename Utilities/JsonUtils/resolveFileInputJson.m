function jsonstruct = resolveFileInputJson(jsonstruct)
    
    fileroot = batmoDir();
    
    if ~isempty(jsonstruct) && isstruct(jsonstruct)
        
        fields_sd = fieldnames(jsonstruct);
        if ismember('isFile', fields_sd)
            filename = jsonstruct.filename;
            fullfilename = fullfile(fileroot, filename);
            jsonsrc = fileread(fullfilename);
            parsedjson = jsondecode(jsonsrc);
            jsonstruct = parsedjson;
        end
    
    end
    
end
