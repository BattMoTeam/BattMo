classdef OxideMembraneControlInputParams < InputParams
    
    properties

        controlType
        
    end
    
    methods
        
        function inputparams = OxideMembraneControlInputParams(jsonstruct)
            
            inputparams = inputparams@InputParams(jsonstruct);
            
        end
        
    end
    
end
