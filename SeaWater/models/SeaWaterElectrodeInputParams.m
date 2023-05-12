classdef SeaWaterElectrodeInputParams < ElectronicComponentInputParams
    
    properties
        
        % Molar Volume
        V
        % structure to describe external coupling
        externalCouplingTerm 
        
    end
    
    methods
        
        function paramobj = SeaWaterElectrodeInputParams(jsonstruct)
            
            paramobj = paramobj@ElectronicComponentInputParams(jsonstruct);
            
        end

    end
    
end
