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
            paramobj.EffectiveElectronicConductivity = 'not used';
        end

    end
    
end
