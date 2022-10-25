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
