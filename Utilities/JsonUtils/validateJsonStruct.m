function is_valid = validateJsonStruct(jsonstruct, reload)

    if nargin == 2 & reload
        reloadModule('validateJsonStruct');
    end

    % Write the json struct to temporary file
    tempfilename = [tempname, '.json'];
    fid = fopen(tempfilename, 'w');
    fprintf(fid, '%s', jsonencode(jsonstruct));
    fclose(fid);

    % Validate
    is_valid = validateJsonFiles({tempfilename});
    
end
