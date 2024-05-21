classdef CO2membraneInputParams < InputParams
    
    properties

        Feed
        Permeate
        permeabilities % structure with permeability
        thickness
        
    end

    methods

        function inputparams = CO2membraneInputParams(jsonstruct)

            inputparams = inputparams@InputParams(jsonstruct);

            inputparams.Feed     = CO2membraneSideInputParams(jsonstruct.Feed);
            inputparams.Permeate = CO2membraneSideInputParams(jsonstruct.Permeate);
            
        end

    end

end
