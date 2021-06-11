classdef ElectrolyteInputParams < ElectroChemicalComponentInputParams

    properties
        
        name
        
        sp
        compnames
        
        sep
        
        conductivityFactor
        
        thermalConductivity
        heatCapacity
        
    end
    
    methods

        function paramobj = ElectrolyteInputParams();
            paramobj = paramobj@ElectroChemicalComponentInputParams();
            paramobj.sp = struct();
            paramobj.compnames = {};
            paramobj.sep = SeparatorInputParams();
            paramobj.EffectiveElectricalConductivity = 'not used';
        end

    end
    
end
