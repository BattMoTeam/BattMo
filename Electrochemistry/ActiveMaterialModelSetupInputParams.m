classdef ActiveMaterialModelSetupInputParams < InputParams
    
    properties
        
        composite
        
        SEImodel
        
    end
    
    methods
        
        function inputparams = ActiveMaterialModelSetupInputParams(jsonstruct)
            
            inputparams = inputparams@InputParams(jsonstruct);
                                          
        end
        
    end
    
    
end

