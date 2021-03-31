classdef ElectroChemicalComponentInputParams < ElectronicComponentInputParams
    
    properties
       
        chargeCarrierName
        chargeCarrierFluxName 
        chargeCarrierSourceName
        chargeCarrierMassConsName
        chargeCarrierAccumName        

    end
    
    methods

        function paramobj = ElectroChemicalComponentInputParams()
            paramobj = paramobj@ElectronicComponentInputParams();
            
            paramobj.chargeCarrierName         = char;
            paramobj.chargeCarrierFluxName     = char;
            paramobj.chargeCarrierSourceName   = char;
            paramobj.chargeCarrierMassConsName = char;
            paramobj.chargeCarrierAccumName    = char;
            
        end
        
    end
    
end
