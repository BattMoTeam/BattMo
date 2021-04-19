classdef ThermoElectronicComponentInputParams < ElectronicComponentInputParams

    properties
        
        thermalConductivity
        heatCapacity
        ohmicResistance
        
    end
    
    methods
        function paramobj = ThermoChemicalComponentInputParams()
            paramobj = paramobj@ElectronicComponentInputParams();
        end
        
    end
    
    
end
