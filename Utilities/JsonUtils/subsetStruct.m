function subjsonstruct = subsetStruct(jsonstruct, fieldnamelists)

    subjsonstruct = [];
    
    for ifs = 1 : numel(fieldnamelists)
        
        fdlist = fieldnamelists{ifs};
        val = getStructField(jsonstruct, fdlist);
        subjsonstruct = setStructField(subjsonstruct, fdlist, val);
        
    end
    
end
