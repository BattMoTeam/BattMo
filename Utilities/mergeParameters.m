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





%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
