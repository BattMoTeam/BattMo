classdef ElectrolyteInputParams < ElectroChemicalComponentInputParams
%
% Input parameter class for :class:`Electrolyte <Electrochemistry.Electrolyte>`
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



%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
