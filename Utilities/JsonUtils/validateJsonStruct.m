function is_valid = validateJsonStruct(jsonstruct, reload)

    if nargin == 2 & reload
        clear classes
        mod = py.importlib.import_module('validationJsonScript');
        py.importlib.reload(mod);
    end

    % Write the json struct to temporary file
    tempfilename = [tempname, '.json'];
    fid = fopen(tempfilename, 'w');
    fprintf(fid, '%s', jsonencode(jsonstruct));
    fclose(fid);

    % Validate
    is_valid = validateJsonFiles({tempfilename});
    
end
