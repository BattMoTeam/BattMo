classdef MaxwellStefanGasDiffusionInputParams < MaxwellStefanDiffusionInputParams
    
    methods

        function inputparams = MaxwellStefanGasDiffusionInputParams(jsonstruct)
            
            inputparams = inputparams@MaxwellStefanDiffusionInputParams(jsonstruct);
            
        end
        
    end

end

