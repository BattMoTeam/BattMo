classdef ThermoElectroChemicalComponentInputParams < ElectroChemicalComponentInputParams
    
    properties
       
        thermalConductivity
        heatCapacity
        ohmicResistance
        
    end
    
    methods

        function paramobj = ThermoElectroChemicalComponentInputParams()
            
            paramobj = paramobj@ElectroChemicalComponentInputParams();
            
        end
        
    end
    
end
