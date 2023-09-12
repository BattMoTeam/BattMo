classdef ProtonicMembraneControlInputParams < InputParams
    
    properties

        controlType
        
    end
    
    methods
        
        function paramobj = ProtonicMembraneControlInputParams(jsonstruct)
            
            paramobj = paramobj@InputParams(jsonstruct);
            
        end
        
    end
    
end
