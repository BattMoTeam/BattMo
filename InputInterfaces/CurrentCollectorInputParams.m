classdef CurrentCollectorInputParams < ElectronicComponentInputParams
    
    properties
        
        volumeFraction 
        electronicConductivity
    end
    
    methods
        
        function params = CurrentCollectorInputParams(params)

            params = params@ElectronicComponentInputParams(params);
            params.EffectiveElectronicConductivity = (params.electronicConductivity) .* (params.volumeFraction).^1.5;
            
        end

    end
    
end
