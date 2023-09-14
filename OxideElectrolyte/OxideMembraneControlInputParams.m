classdef OxideMembraneControlInputParams < InputParams
    
    properties

        controlType
        
    end
    
    methods
        
        function paramobj = OxideMembraneControlInputParams(jsonstruct)
            
            paramobj = paramobj@InputParams(jsonstruct);
            
        end
        
    end
    
end
