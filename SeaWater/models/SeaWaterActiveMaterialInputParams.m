classdef SeaWaterActiveMaterialInputParams < ComponentInputParams

    properties
        
        chargeCarrierName
        stochElectron
        etaMax
    end
    
    methods

        function paramobj = SeaWaterActiveMaterialInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
        end
        
    end
    

end
