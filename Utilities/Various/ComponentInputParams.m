classdef ComponentInputParams < InputParams
    properties
        G % grid
    end
    
    methods

        function paramobj = ComponentInputParams(jsonstruct)
            
            paramobj = paramobj@InputParams(jsonstruct);
            
        end
        
    end
    
end
