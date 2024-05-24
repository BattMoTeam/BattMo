classdef ActiveMaterialModelSetupInputParams < InputParams
    
    properties
        
        composite
        
        SEImodel
        
    end
    
    methods
        
        function inputparams = ActiveMaterialModelSetupInputParams(jsonstruct)

            inputparams = inputparams@InputParams(jsonstruct);
            
            inputparams = inputparams.setupDefault();
        
        end
        
        function inputparams = setupDefault(inputparams)
        
            inputparams = setupDefault@InputParams(inputparams);
            
            inputparams = inputparams.setDefault({'composite'}, false);
            inputparams = inputparams.setDefault({'SEImodel'}, 'none');
            
        end

    end
    
    
end

