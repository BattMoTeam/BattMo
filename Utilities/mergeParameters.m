function paramobj = mergeParameters(paramobj, paramnames, varargin)

% paramnames is a list of parameter names. The function assign all empty values to the given non-empty ones

% If all parameter valus are empty, nothing is done
    
% If 'force' is set to true, the first non-empty parameter is used to enforce consistancy
    
    opt = struct('force', false); 
    opt = merge_options(opt, varargin{:});
    
    vals = cellfun(@(paramname) paramobj.getParam(paramname), paramnames, 'uniformoutput', false);

    ind = cellfun(@(val) ~isempty(val), vals);

    if any(nnz(ind))

        iref = find(ind);
        iref = iref(1);
        val = vals{iref};
        
        if ~opt.force & (nnz(ind) > 1)
            % We check that all the given values are consistent
            vals = vals(ind);
            if ischar(val)
                vals = unique(vals);
            elseif islogical(val)| isnumeric(val)
                vals = unique([vals{:}]);
            else
                error('paramter type is not supported by this function');
            end
            assert(numel(vals) == 1, 'inconsistant parameters have been given')
                
        end

        for iparam = 1 : numel(paramnames)
            paramname = paramnames{iparam};
            paramobj = paramobj.setParam(paramname, val);
        end
        
    end

    
end


