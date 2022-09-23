function paramobj = mergeParameters(paramobj, paramname1, paramname2, varargin)

    opt = struct('force', true); 
    opt = merge_options(opt, varargin{:});
    
    val1 = paramobj.getParam(paramname1);
    val2 = paramobj.getParam(paramname2);

    if isempty(val1) && ~isempty(val2)
        paramobj = paramobj.setParam(paramname1, val2);
    elseif ~isempty(val1) && isempty(val2)
        paramobj = paramobj.setParam(paramname2, val1);
    elseif ~isempty(val1) && ~isempty(val2)
        if opt.force
            % we use value given from the first one.
            paramobj = paramobj.setParam(paramname2, val1);
        else
            assert(val1 == val2, 'inputs are not consistent');
        end
    end
    
end


