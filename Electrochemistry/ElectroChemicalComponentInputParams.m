classdef ElectroChemicalComponentInputParams < ElectronicComponentInputParams
%
% Input class for :class:`ElectronicComponent <Electrochemistry.ElectronicComponent>`
%
    properties
       
        chargeCarrierName

    end
    
    methods

        function paramobj = ElectroChemicalComponentInputParams(jsonstruct)
            paramobj = paramobj@ElectronicComponentInputParams(jsonstruct);
        end
        
    end
    
end
