function params = configToModelParameter(setup, config)

    % Remove keys not for ModelParameter
    keys = config.Properties.VariableNames;
    [keys, keysIdx] = filter(keys);
    n = numel(keys);

    opts = cell(2*n, 1);
    opts(1:2:2*n) = keys;

    numVars = numel(config.Row);
    params = cell(numVars, 1);

    for k = 1:numVars
        % Remove keys
        conf = config(k, :);
        vals = table2cell(conf);
        vals = filter(vals, keysIdx);

        % Save
        opts(2:2:2*n) = vals;

        % Store parameters
        param = ModelParameter(setup, opts{:});
        params{k} = param;
    end

end

function [b, ia] = filter(a, idx)

    if nargin == 1
        props = properties('ModelParameter');
        [b, ia] = intersect(a, props);
    else
        b = a(idx);
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
