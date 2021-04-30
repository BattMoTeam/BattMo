classdef ThermoElectroChemicalComponentInputParams < ElectroChemicalComponentInputParams
    
    properties
       
        thermalConductivity
        heatCapacity
        
    end
    
    methods

        function paramobj = ThermoElectroChemicalComponentInputParams()
            
            paramobj = paramobj@ElectroChemicalComponentInputParams();
            
        end
        
    end
    
end
