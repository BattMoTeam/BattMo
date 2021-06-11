classdef LithiumElectrodeInputParams < ElectrodeInputParams
    
    methods

        function paramobj = LithiumElectrodeInputParams()
            
            paramobj = paramobj@ElectrodeInputParams();
            paramobj.eac.chargeCarrierName         = 'Li';
            
        end

    end
    
end
