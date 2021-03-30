classdef ElectrolyteInputParams < ActiveElectroChemicalComponentInputParams
    
    properties
        electrolyteType
    end
    
    methods
        
        function paramobj = setup(paramobj, params)
        % params struct should contain valid fields for ActiveElectroChemicalComponentInputParams
        %
        % and fields
        %
        % - electrolyteType
        %
            
            paramobj = setup@ActiveElectroChemicalComponent(paramobj, params);
            
        end
        
    end
    
end
