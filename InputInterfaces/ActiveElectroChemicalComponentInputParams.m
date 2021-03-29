classdef ActiveElectroChemicalComponentInputParams < ElectroChemicalComponentInputParams
% for the moment, not much here
    
    methods
        
        function paramobj = ActiveElectroChemicalComponentInputParams(params)
        % params struct should contain valid fields for ElectroChemicalComponentInputParams
            paramobj = paramobj@ElectroChemicalComponentInputParams(params);
        end
        
    end
    
end
