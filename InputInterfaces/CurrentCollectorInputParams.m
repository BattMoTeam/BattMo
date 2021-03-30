classdef CurrentCollectorInputParams < ElectronicComponentInputParams
    
    properties
        volumeFraction 
        electronicConductivity
    end
    
    methods

        function paramobj = setup(paramobj, params)
        %
        % params should contain fields for ComponentInputParams and 
        % - volumeFraction
        % - electronicConductivity
            
            
            volumeFraction = params.volumeFraction;
            electronicConductivity = params.electronicConductivity;
            params.EffectiveElectronicConductivity = (electronicConductivity) .* (volumeFraction).^1.5;
            
            paramobj = setup@ElectronicComponentInputParams(paramobj, params);

            paramobj.volumeFraction = volumeFraction;
            paramobj.electronicConductivity = electronicConductivity;
            
        end
        
    end
    
end
