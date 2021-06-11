classdef ElectroChemicalComponentInputParams < ElectronicComponentInputParams
    
    properties
       
        chargeCarrierName

    end
    
    methods

        function paramobj = ElectroChemicalComponentInputParams()

            paramobj = paramobj@ElectronicComponentInputParams();
            paramobj.chargeCarrierName = char();
            
        end
        
    end
    
end
