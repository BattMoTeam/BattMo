classdef ThermoElectronicComponentInputParams < ElectronicComponentInputParams

    properties
        
        thermalConductivity
        heatCapacity % in [J][K]^-1[m]^-3
        
    end
    
    methods
        function paramobj = ThermoElectronicComponentInputParams()
            paramobj = paramobj@ElectronicComponentInputParams();
        end
        
    end
    
    
end
