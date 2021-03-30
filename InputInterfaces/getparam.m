function val = getparam(params, fdname, default)
    if nargin > 2
        default = [];
    end
    
    if isfield(params, fd)
        val = params.(fdname)
    else
        val = default;
    end
    
end

