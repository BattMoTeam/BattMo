classdef ElectrodeInputParams < ComponentInputParams
%
% Input parameter class for :class:`Electrode <Electrochemistry.Electrodes.Electrode>`
%        
    properties
        
        %
        % Input parameter for the electrode active component (:class:`ActiveMaterialInputParams
        % <Electrochemistry.ActiveMaterialInputParams>`)
        %
        ActiveMaterial
        
        
        %
        % Input parameter for the current collector (:class:`CurrentCollectorInputParams
        % <Electrochemistry.CurrentCollectorInputParams>`)
        %
        CurrentCollector
        
        %
        % Coupling term specification
        %
        couplingTerm
        
    end
    
    methods

        function paramobj = ElectrodeInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);

            am = 'ActiveMaterial';
            cc  = 'CurrentCollector';
            
            pick = @(fd) pickField(jsonstruct, fd);
            paramobj.(am) = ActiveMaterialInputParams(pick(am));
            paramobj.(cc) = CurrentCollectorInputParams(pick(cc));
            
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
