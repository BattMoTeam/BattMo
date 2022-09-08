function paramobj = mergeParameters(paramobj, paramname1, paramname2)

    val1 = paramobj.getParam(paramname1);
    val2 = paramobj.getParam(paramname2);

    if isempty(val1) && ~isempty(val2)
        paramobj = paramobj.setParam(paramname1, val2);
    elseif ~isempty(val1) && isempty(val2)
        paramobj = paramobj.setParam(paramname1, val1);
    elseif ~isempty(val1) && ~isempty(val2)
        %% TODO : add more precise error message that indicates which parameters are not consistent
        assert(val1 == val2, 'inputs are not consistent')
    end
    
end


