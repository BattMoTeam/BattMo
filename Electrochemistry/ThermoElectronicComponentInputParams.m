classdef ThermoElectronicComponentInputParams < ElectronicComponentInputParams

    properties
        
        thermalConductivity
        heatCapacity
        
    end
    
    methods
        function paramobj = ThermoElectronicComponentInputParams()
            paramobj = paramobj@ElectronicComponentInputParams();
        end
        
    end
    
    
end
