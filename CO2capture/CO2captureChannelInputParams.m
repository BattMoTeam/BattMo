classdef CO2captureChannelInputParams < ComponentInputParams

    properties


        couplingTerms
        
        Control
        
        gasSpecies

        rateCoefficient

    end

    methods
        
        function inputparams = CO2captureChannelInputParams(jsonstruct)
            
            inputparams = inputparams@ComponentInputParams(jsonstruct);

        end
        
    end

end
        
