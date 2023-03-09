function is_valid = validateJsonStruct(jsonstruct)

    % Write the json struct to a temporary file
    tempfilename = [tempname, '.json'];
    fid = fopen(tempfilename, 'w');
    fprintf(fid, '%s', jsonencode(jsonstruct));
    fclose(fid);

    % Validate
    is_valid = validateJsonFiles({tempfilename});
    
end
