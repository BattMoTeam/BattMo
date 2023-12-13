classdef MagnesiumActiveMaterialInputParams < SeaWaterActiveMaterialInputParams

    properties

        Asurf

    end
    
    methods

        function inputparams = MagnesiumActiveMaterialInputParams(jsonstruct)
    
            inputparams = inputparams@SeaWaterActiveMaterialInputParams(jsonstruct);
            
        end
        
    end
    
end

