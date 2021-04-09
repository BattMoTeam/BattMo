classdef LithiumElectrodeInputParams < ElectrodeInputParams
    
    methods

        function paramobj = LithiumElectrodeInputParams()
            
            paramobj = paramobj@ElectrodeInputParams();

            paramobj.eac.chargeCarrierName         = 'Li';
            paramobj.eac.chargeCarrierFluxName     = 'LiFlux';
            paramobj.eac.chargeCarrierSourceName   = 'LiSource';
            paramobj.eac.chargeCarrierMassConsName = 'massCons';
            paramobj.eac.chargeCarrierAccumName    = 'LiAccum';
            
        end

    end
    
end
