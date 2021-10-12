function paramobj = assignJsonParams(paramobj, jsonstruct)
    
    jsonstruct = resolveFileInput(jsonstruct);
    fields_pobj = properties(paramobj);

    for ind = 1 : numel(fields_pobj)
        
        fd = fields_pobj{ind};
        if isfield(jsonstruct, fd)
            paramobj.(fd) = jsonstruct.(fd);
        end
        
    end
    
end

function jsonstruct = resolveFileInput(jsonstruct)
    
    fileroot = batmoDir();
    
    if ~isempty(jsonstruct)
        
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
