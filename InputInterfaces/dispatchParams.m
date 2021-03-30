function paramobj = dispatchParams(paramobj, params, fdnames)
    if iscell(fdnames)
        for ind = 1 : numel(fdnames)
            fdname = fdnames{ind};
            paramobj.(fdname) = getparam(params, fdname);
        end
    elseif ischar(fdnames)
        paramobj.(fdnames) = getparam(params, fdnames);        
    else
        error('type of fdnames not recognized');
    end
        
end
