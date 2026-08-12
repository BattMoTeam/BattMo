classdef CO2captureInputParams < InputParams
    
    properties

        Feed
        Permeate
        permeabilities % structure with permeability
        thickness
        
    end

    methods

        function inputparams = CO2captureInputParams(jsonstruct)

            inputparams = inputparams@InputParams(jsonstruct);

            inputparams.Feed     = CO2captureChannelInputParams(jsonstruct.Feed);
            inputparams.Permeate = CO2captureChannelInputParams(jsonstruct.Permeate);
            
        end

    end

end
