function result = isEqualStructField(jsonstruct, fieldnamelist, testValue)

    if isUnAssigned(jsonstruct, fieldnamelist)
        result = false;
        return
    end
    
    value = getStructField(jsonstruct, fieldnamelist);
    result = isequal(value, testValue);
    
end
