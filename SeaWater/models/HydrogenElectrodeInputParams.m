classdef HydrogenElectrodeInputParams < SeaWaterElectrodeInputParams
    
    methods
        
        function paramobj = HydrogenElectrodeInputParams(jsonstruct)
            
            paramobj = paramobj@SeaWaterElectrodeInputParams(jsonstruct);
            
        end

    end
    
end
