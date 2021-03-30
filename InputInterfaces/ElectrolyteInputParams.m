classdef ElectrolyteInputParams < ElectroChemicalComponentInputParams

    properties
        
        sp
        compnames
        ncomp
        
    end
    
    methods
        
        function paramobj = setup(paramobj, params)
            paramobj = setup@ElectroChemicalComponent(paramobj, params);
        end
        
    end
    
end
