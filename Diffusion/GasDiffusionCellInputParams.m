classdef GasDiffusionCellInputParams < MaxwellStefanGasDiffusionInputParams

    properties

        Control
        
    end
    
    methods

        function inputparams = GasDiffusionCellInputParams(jsonstruct)

            jsonstruct = equalizeJsonStructFields(jsonstruct, {'compNames', {'Control', 'compNames'}});
            
            inputparams = inputparams@MaxwellStefanGasDiffusionInputParams(jsonstruct);

            inputparams.Control = GasDiffusionCellControlInputParams(jsonstruct.Control);
            
        end
        
    end

end

