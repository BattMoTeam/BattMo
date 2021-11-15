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
        updateDiffusionCoefficientFunc
        
        density % [kg m^-3] (Note : only of the liquid part, the density of the separator is given there)
        
    end
    
    methods

        function paramobj = ElectrolyteInputParams(jsonstruct);
            
            paramobj = paramobj@ElectroChemicalComponentInputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);
            paramobj.Separator = SeparatorInputParams(pick('Separator'));
            paramobj.updateConductivityFunc = FunctionInputParams(pick('updateConductivityFunc'));
            paramobj.updateDiffusionCoefficientFunc = FunctionInputParams(pick('updateDiffusionCoefficientFunc'));
            
            paramobj.EffectiveElectricalConductivity = 'not used';
        end

    end
    
end
