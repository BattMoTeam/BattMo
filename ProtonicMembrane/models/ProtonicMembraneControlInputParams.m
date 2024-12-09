classdef ProtonicMembraneControlInputParams < InputParams
    
    properties

        controlType
        Imax
        
    end
    
    methods
        
        function inputparams = ProtonicMembraneControlInputParams(jsonstruct)

            jsonstruct = setDefaultJsonStructField(jsonstruct, {'Imax'}, 0);
            inputparams = inputparams@InputParams(jsonstruct);
            
        end
        
    end
    
end
