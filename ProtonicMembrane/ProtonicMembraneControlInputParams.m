classdef ProtonicMembraneControlInputParams < InputParams
    
    properties

        controlType
        
    end
    
    methods
        
        function inputparams = ProtonicMembraneControlInputParams(jsonstruct)
            
            inputparams = inputparams@InputParams(jsonstruct);
            
        end
        
    end
    
end
