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
        
        electrolyteType
        
        updateConductivityFunc
        
    end
    
    methods

        function paramobj = ElectrolyteInputParams(jsonstruct);
            
            paramobj = paramobj@ElectroChemicalComponentInputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);
            paramobj.Separator = SeparatorInputParams(pick('Separator'));
            paramobj.updateConductivityFunc = FunctionInputParams(pick('updateConductivityFunc'));
            
            paramobj.EffectiveElectricalConductivity = 'not used';
        end

    end
    
end
