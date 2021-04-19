classdef ThermoElectroChemicalComponentInputParams < ElectroChemicalComponentInputParams
    
    properties
       
        ThermalConductivity

    end
    
    methods

        function paramobj = ThermoElectroChemicalComponentInputParams()
            
            paramobj = paramobj@ElectroChemicalComponentInputParams();
            
            
        end
        
    end
    
end
