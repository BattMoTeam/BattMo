classdef MagnesiumActiveMaterialInputParams < SeaWaterActiveMaterialInputParams

    properties

        Asurf

    end
    
    methods

        function paramobj = MagnesiumActiveMaterialInputParams(jsonstruct)
    
            paramobj = paramobj@SeaWaterActiveMaterialInputParams(jsonstruct);
            
        end
        
    end
    
end

