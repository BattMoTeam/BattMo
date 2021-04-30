classdef ThermoElectroChemicalComponentInputParams < ElectroChemicalComponentInputParams
    
    properties
       
        thermalConductivity
        heatCapacity % in [J][K]^-1[m]^-3
        
    end
    
    methods

        function paramobj = ThermoElectroChemicalComponentInputParams()
            
            paramobj = paramobj@ElectroChemicalComponentInputParams();
            
        end
        
    end
    
end
