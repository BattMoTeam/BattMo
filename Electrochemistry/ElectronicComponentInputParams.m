classdef ElectronicComponentInputParams < ComponentInputParams
%
% Input class for :class:`ElectronicComponent <Electrochemistry.ElectronicComponent>`
%    
    properties
        EffectiveElectricalConductivity;
    end
    
    
    methods

        function paramobj = ElectronicComponentInputParams(jsonstruct)
            paramobj = paramobj@ComponentInputParams(jsonstruct);
        end
        
    end
    

    
end
