classdef ElectrolyteInputParams < ElectroChemicalComponentInputParams
%
% Input class for :class:`Electrolyte <Electrochemistry.Electrolyte>`
%    
    properties
        
        name
        
        sp
        compnames
        
        Separator
        
        conductivityFactor
        
        thermalConductivity
        heatCapacity
        
    end
    
    methods

        function paramobj = ElectrolyteInputParams();
            paramobj = paramobj@ElectroChemicalComponentInputParams();
            paramobj.sp = struct();
            paramobj.compnames = {};
            paramobj.Separator = SeparatorInputParams();
            paramobj.EffectiveElectricalConductivity = 'not used';
        end

    end
    
end
