classdef GasDiffusionCellInputParams < MaxwellStefanGasDiffusionInputParams

    properties

        Control

        externalCouplingTerms % Cell array with external coupling terms
        
    end
    
    methods

        function inputparams = GasDiffusionCellInputParams(jsonstruct)

            jsonstruct = equalizeJsonStructFields(jsonstruct, {'componentNames', {'Control', 'componentNames'}});
            
            inputparams = inputparams@MaxwellStefanGasDiffusionInputParams(jsonstruct);

            inputparams.Control = GasDiffusionCellControlInputParams(jsonstruct.Control);
            
        end
        
    end

end

