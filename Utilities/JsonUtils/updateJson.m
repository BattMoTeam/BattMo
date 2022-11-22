function json = updateJson(json, params, varargin)

    opt = struct('validate', true, ...
                 'tempfilename', []);
    opt = merge_options(opt, varargin{:});

    assert(size(params,1)==1 || size(params,2)==1, 'params should be a cell array with one column or one row');
    assert(rem(numel(params), 2) == 0, 'params should be contain key and value pairs');

    % Extract
    keys = params(1:2:end);
    vals = params(2:2:end);

    % Update struct
    for k = 1:numel(keys)
        key = split(keys{k}, '.');
        val = vals{k};
        json = setfield(json, key{:}, val);
    end

    if opt.validate
        % Write the json struct to file
        if isempty(opt.tempfilename)
            tempfilename = [tempname, '.json'];
        else
            tempfilename = opt.tempfilename;
        end
        fid = fopen(tempfilename, 'w');
        fprintf(fid, '%s', jsonencode(json));
        fclose(fid);

        % Validate
        is_valid = py.validationJsonScript.validate(tempfilename);
        assert(is_valid);
    end

end
