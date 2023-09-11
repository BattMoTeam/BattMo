classdef ProtonicMembraneControlInputParams < InputParams
    
    properties

        controlType
        controlValue
        
    end
    
    methods
        
        function paramobj = ProtonicMembraneControlInputParams(jsonstruct)
            
            paramobj = paramobj@InputParams(jsonstruct);
            
        end
        
    end
    
end
