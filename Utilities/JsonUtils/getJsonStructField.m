function value = getJsonStructField(jsonstruct, fieldnamelist)

    if ischar(fieldnamelist)
        % handle case where fieldnamelist is just a char
        fieldnamelist = {fieldnamelist};
        value = getJsonStructField(jsonstruct, fieldnamelist);
        return
    end

    fieldname = fieldnamelist{1};

    if numel(fieldnamelist) > 1

        fieldnamelist = fieldnamelist(2 : end);
        getValue = false;
        
    else
        
        getValue = true;
        
    end

    if ~isfield(jsonstruct, fieldname)

        value = UnAssigned();

    else

        if getValue

            value = jsonstruct.(fieldname);

        else

            value = getJsonStructField(jsonstruct.(fieldname), fieldnamelist);

        end
        
    end
        
    
end
