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

        function paramobj = SeaWaterElectrolyteInputParams(jsonstruct)
            
            paramobj = paramobj@SeaWaterElectrolyteNoPrecipitationInputParams(jsonstruct);
            
        end
        
    end
end
