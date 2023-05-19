classdef HydrogenActiveMaterialInputParams < SeaWaterActiveMaterialInputParams

    properties
        
        Asurf % 
        
    end

    methods

        function paramobj = HydrogenActiveMaterialInputParams(jsonstruct)
            
            paramobj = paramobj@SeaWaterActiveMaterialInputParams(jsonstruct);
            
        end
        
    end
    
end
