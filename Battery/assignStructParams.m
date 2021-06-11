function paramobj = assignStructParams(paramobj, structdata)
    
    fields_sd = fieldnames(structdata);
    fields_pobj = fieldnames(paramobj);
    
    for ind = 1 : numel(fields_sd)
        
        fd = fields_sd{ind};
        
        % if isclass(paramobj)
            % if paramobj is a class, we check here that it the field fd matches a property of the class
            % assert(ismember(fd, fields_pobj), 'field in input data is not recognized');
        % end
        
        if isnumeric(structdata.(fd))
            paramobj.(fd) = structdata.(fd);
        elseif ischar(structdata.(fd))
            paramobj.(fd) = structdata.(fd);
        elseif iscell(structdata.(fd))
            paramobj.(fd) = structdata.(fd);
        elseif isstruct(structdata.(fd))
            paramobj.(fd) = assignStructParams(paramobj.(fd), structdata.(fd));
        end

    end
    
end

