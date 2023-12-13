classdef SeaWaterElectrodeInputParams < ElectronicComponentInputParams
    
    properties
        
        % Molar Volume
        V
        % structure to describe external coupling
        externalCouplingTerm 
        
    end
    
    methods
        
        function inputparams = SeaWaterElectrodeInputParams(jsonstruct)
            
            inputparams = inputparams@ElectronicComponentInputParams(jsonstruct);
            
        end

    end
    
end
