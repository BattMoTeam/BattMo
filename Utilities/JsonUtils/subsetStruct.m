function subjsonstruct = subsetJsonStruct(jsonstruct, fieldnamelists)

    subjsonstruct = [];
    
    for ifs = 1 : numel(fieldnamelists)
        
        fdlist = fieldnamelists{ifs};
        val = getJsonStructField(jsonstruct, fdlist);
        subjsonstruct = setJsonStructField(subjsonstruct, fdlist, val);
        
    end
    
end
