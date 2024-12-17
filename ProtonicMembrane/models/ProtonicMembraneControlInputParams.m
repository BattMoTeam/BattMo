classdef ProtonicMembraneControlInputParams < InputParams
    
    properties

        controlType
        I
        
    end
    
    methods
        
        function inputparams = ProtonicMembraneControlInputParams(jsonstruct)

            jsonstruct = setDefaultJsonStructField(jsonstruct, {'I'}, 0);
            inputparams = inputparams@InputParams(jsonstruct);
            
        end
        
    end
    
end
