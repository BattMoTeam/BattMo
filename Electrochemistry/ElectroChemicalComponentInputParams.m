classdef ElectroChemicalComponentInputParams < ElectronicComponentInputParams
%
% Input class for :class:`ElectronicComponent <Electrochemistry.ElectronicComponent>`
%
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
