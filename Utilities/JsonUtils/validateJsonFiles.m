function is_valid = validateJsonFiles(jsonfiles, reload)

    if ~iscell(jsonfiles)
        jsonfiles = {jsonfiles};
    end

    if nargin == 2 & reload
        reloadModule('validationJsonScript');
    end

    % Validate using python script
    for k = 1:numel(jsonfiles)
        jsonfile = jsonfiles{k};
        is_valid{k} = py.validationJsonScript.validate(jsonfile);
        assert(is_valid{k}, 'jsonfile %s is not valid', jsonfile);
    end

    if numel(jsonfiles) == 1
        is_valid = is_valid{1};
    end

end
