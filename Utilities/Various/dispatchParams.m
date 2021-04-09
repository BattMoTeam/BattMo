function obj = dispatchParams(obj, params, fdnames)
    if iscell(fdnames)
        for ind = 1 : numel(fdnames)
            fdname = fdnames{ind};
            obj = dispatchParams(obj, params, fdname);
        end
    elseif ischar(fdnames)
        fdname = fdnames;
        if isfield(params, fdname) | isprop(params, fdname)
            obj.(fdname) = params.(fdname);
        end
    else
        error('type of fdnames not recognized');
    end
        
end
