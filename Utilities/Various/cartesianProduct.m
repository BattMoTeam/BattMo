function [cartinds, cartvalues] = cartesianProduct(varargin)

%
%
% Each argument in varargin is a (cell-)array of values (any type such as string, boolean or scalar).
%
% The function setup a *cartesian products* of all these values. It means that the output cartinds is an array of dimension
% (number of combination)x(number of cells in params) such as, for each combination given by a row i, cartind(i, j) gives the
% index in params{j} that corresponds to this combination.
%
% If the vargout cartvalues is given, the function also sets a cell-array cartvalues where cartvalues{i} is a structure field names
% given by the varargin name (we use function inputname). For a combination given by the index i, cartvalues{i} gives the
% *values* of the (cell-)array sent in vararing
%     
% The following example should help to understand the explanations above. Run and inspect the output
%
%    a = {true, false};
%    b = [1; 2; 3];
%    c = {'a', 'b'};
%    
%    [cartind, cartvalues] = cartesianProduct(a, b, c);
    

    genInputs = (nargout > 1);

    for iarg = 1 : nargin
        inds{iarg} = 1 : numel(varargin{iarg});
        if genInputs
            names{iarg} = inputname(iarg);
        end
    end

    [inds{:}] = ndgrid(inds{:});

    for iinds = 1 : numel(inds)
        
        cartinds(:, iinds) = reshape(inds{iinds}, [], 1);

    end

    if genInputs

        for icomb = 1 : size(cartinds, 1)
            
            clear cartvalue
            for iarg = 1 : nargin
                ind = cartinds(icomb, iarg);
                if iscell(varargin{iarg})
                    val = varargin{iarg}{ind};
                else
                    val = varargin{iarg}(ind);
                end
                cartvalue.(names{iarg}) = val;
            end
            cartvalues{icomb} = cartvalue;
        end
        
    end
    
end
