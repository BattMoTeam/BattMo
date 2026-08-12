classdef CO2captureChannelInputParams < ComponentInputParams

    properties

        Control
        
        gasSpecies

        rateCoefficient
        
        area
        length
        

    end

    methods
        
        function inputparams = CO2captureChannelInputParams(jsonstruct)
            
            inputparams = inputparams@ComponentInputParams(jsonstruct);

        end
        
    end

end
        
