function paramobj = assignJsonParams(paramobj, jsonstruct)
    
    jsonstruct = resolveFileInputJson(jsonstruct);
    fields_pobj = properties(paramobj);

    for ind = 1 : numel(fields_pobj)
        
        fd = fields_pobj{ind};
        if isfield(jsonstruct, fd)
            paramobj.(fd) = jsonstruct.(fd);
        end
        
    end
    
end
