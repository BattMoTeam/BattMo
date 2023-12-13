classdef SeaWaterActiveMaterialInputParams < ComponentInputParams

    properties
        
        chargeCarrierName
        stochElectron
        etaMax
    end
    
    methods

        function inputparams = SeaWaterActiveMaterialInputParams(jsonstruct)

            inputparams = inputparams@ComponentInputParams(jsonstruct);
            
        end
        
    end
    

end
