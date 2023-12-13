classdef HydrogenElectrodeInputParams < SeaWaterElectrodeInputParams
    
    methods
        
        function inputparams = HydrogenElectrodeInputParams(jsonstruct)
            
            inputparams = inputparams@SeaWaterElectrodeInputParams(jsonstruct);
            
        end

    end
    
end
