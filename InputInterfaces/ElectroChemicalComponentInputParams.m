classdef ElectroChemicalComponentInputParams < ElectronicComponentInputParams
    
    properties
       
        chargeCarrierName
        chargeCarrierFluxName 
        chargeCarrierSourceName
        chargeCarrierMassConsName
        chargeCarrierAccumName        

    end
    
    methods
        
        function paramobj = setup(parmobj, params)
        % params struct should contain valid fields for ElectronicComponentInputParams and the fields
        % - chargeCarrierName
        % - chargeCarrierFluxName 
        % - chargeCarrierSourceName
        % - chargeCarrierMassConsName
        % - chargeCarrierAccumName        
            
            fdnames = {'chargeCarrierName', ...
                       'chargeCarrierFluxName', ...
                       'chargeCarrierSourceName', ...
                       'chargeCarrierMassConsName', ...
                       'chargeCarrierAccumName'};
            paramobj = dispatchParams(paramobj, params, fdnames);

            paramobj = setup@ElectronicComponentInputParams(paramobj, params);
            
        end
        
    end
    
end
