classdef ElectrolyteInputParams < ActiveElectroChemicalComponentInputParams
    
    properties
        electrolyteType
    methods
        
        function paramobj = ElectrolyteInputParams(params)
        % params struct should contain valid fields for ActiveElectroChemicalComponentInputParams
        %
        % and fields
        %
        % - electrolyteType
        %
            
            paramobj = paramobj@ActiveElectroChemicalComponent(params);
            
        end
        
        
        
        
    end
    
end
