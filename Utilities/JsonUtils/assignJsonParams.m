function paramobj = assignJsonParams(paramobj, jsonstruct)
    
    paramobjFds = properties(paramobj);

    for ind = 1 : numel(paramobjFds)
        
        fd = paramobjFds{ind};
        if isfield(jsonstruct, fd)
            paramobj.(fd) = jsonstruct.(fd);
        end
        
    end
    
end
