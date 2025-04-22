function result = isEqualJsonStructField(jsonstruct, fieldnamelist, testValue)

    if isUnAssigned(jsonstruct, fieldnamelist)
        result = false;
        return
    end
    
    value = getJsonStructField(jsonstruct, fieldnamelist);
    result = isequal(value, testValue);
    
end
