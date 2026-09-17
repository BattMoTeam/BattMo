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
    end

    [inds{:}] = ndgrid(inds{:});

    for iinds = 1 : numel(inds)
        cartinds(:, iinds) = reshape(inds{iinds}, [], 1);
    end

    if genInputs

        for iarg = 1 : nargin
            names{iarg} = inputname(iarg);
        end

        cartvalues = cell(size(cartinds, 1), 1);

        for icomb = 1 : size(cartinds, 1)

            for iarg = 1 : nargin
                ind = cartinds(icomb, iarg);
                if iscell(varargin{iarg})
                    val = varargin{iarg}{ind};
                else
                    val = varargin{iarg}(ind);
                end
                cartvalues{icomb}.(names{iarg}) = val;
            end

        end

    end

end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
