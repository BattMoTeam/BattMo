classdef ActiveElectroChemicalComponentInputParams < ElectroChemicalComponentInputParams
% for the moment, not much here
    
    methods
        
        function params = ActiveElectroChemicalComponentInputParams(params)
            params = params@ElectroChemicalComponentInputParams(params);
        end
        
    end
    
end
