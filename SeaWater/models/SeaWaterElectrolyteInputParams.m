classdef SeaWaterElectrolyteInputParams < SeaWaterElectrolyteNoPrecipitationInputParams

    properties

        dischargeProductMolarVolume
        superOversaturationRatio
        solidPrecipitatePorosity
        characteristicPoreDiameter
        nucleationMaximum        
        nucleationRate
        nucleationActivation
        
    end

    methods

        function inputparams = SeaWaterElectrolyteInputParams(jsonstruct)
            
            inputparams = inputparams@SeaWaterElectrolyteNoPrecipitationInputParams(jsonstruct);
            
        end
        
    end
end
