classdef ProtonicMembraneControlInputParams < InputParams
    
    properties

        controlType
        useCurrentDensity % If true, then we use the value of currentDensity. Otherwise, a value for I the total current should
                          % be given.
        currentDensity
        I
        
        area % if current density is used, we need the area. This value will typically be assigned by the grid constructor
        
    end
    
    methods
        
        function inputparams = ProtonicMembraneControlInputParams(jsonstruct)

            jsonstruct  = setDefaultJsonStructField(jsonstruct, {'useCurrentDensity'}, false);
            jsonstruct  = setDefaultJsonStructField(jsonstruct, {'I'}, 0);
            inputparams = inputparams@InputParams(jsonstruct);
            
        end
        
    end
    
end
