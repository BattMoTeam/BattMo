classdef ThermoElectronicComponentInputParams < ElectronicComponentInputParams

    properties
        
        thermalConductivity
        heatCapacity
        ohmicResistance
        
    end
    
    methods
        function paramobj = ThermoElectronicComponentInputParams()
            paramobj = paramobj@ElectronicComponentInputParams();
        end
        
    end
    
    
end
