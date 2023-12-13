classdef HydrogenActiveMaterialInputParams < SeaWaterActiveMaterialInputParams

    properties
        
        Asurf % 
        
    end

    methods

        function inputparams = HydrogenActiveMaterialInputParams(jsonstruct)
            
            inputparams = inputparams@SeaWaterActiveMaterialInputParams(jsonstruct);
            
        end
        
    end
    
end
