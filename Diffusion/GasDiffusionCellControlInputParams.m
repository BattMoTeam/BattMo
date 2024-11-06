classdef GasDiffusionCellControlInputParams < InputParams

    properties
        
        % Cell array with the component names
        compNames

        
    end
    
    methods

        function inputparams = GasDiffusionCellControlInputParams(jsonstruct)

            inputparams = inputparams@InputParams(jsonstruct);
            
        end

        
    end
    

end
