classdef ElectrolyteInputParams < ElectroChemicalComponentInputParams

    properties
        
        sp
        compnames
        ncomp
        
        sep
        
    end
    
    methods

        function paramobj = ElectrolyteInputParams();
            paramobj = paramobj@ElectroChemicalComponentInputParams();
            paramobj.sp = struct();
            paramobj.compnames = {};
            paramobj.sep = SeparatorInputParams();
        end
        
        function paramobj = setup(paramobj, params)
            paramobj = setup@ElectroChemicalComponent(paramobj, params);
        end
        
    end
    
end
