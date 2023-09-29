classdef ElectrodeInputParams < ComponentInputParams
%
% Input parameter class for :code:`Electrode` model
%        
    properties
        
        %% Sub-Models
        
        Coating
        CurrentCollector
        
        %% Coupling term specification
        
        couplingTerm
        
        %% Parameters assigned at setup
        
        include_current_collectors
        use_thermal
        
    end
    
    methods

        function paramobj = ElectrodeInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);

            co = 'Coating';
            cc = 'CurrentCollector';
            
            pick = @(fd) pickField(jsonstruct, fd);

            paramobj.(co) = CoatingInputParams(pick(co));
            paramobj.(cc) = CurrentCollectorInputParams(pick(cc));

            if isempty(paramobj.(cc))
                paramobj.include_current_collectors = false;
            else
                paramobj.include_current_collectors = true;
            end
            
            paramobj = paramobj.validateInputParams();

        end

        function paramobj = validateInputParams(paramobj)

            am = 'Coating';
            cc ='CurrentCollector';
            
            paramobj = mergeParameters(paramobj, {{'use_thermal'}, {am, 'use_thermal'}});
            paramobj.(am) = paramobj.(am).validateInputParams();

            if ~isempty(paramobj.include_current_collectors)
                paramobj = mergeParameters(paramobj, {{'use_thermal'}, {cc, 'use_thermal'}});
                paramobj.(cc) = paramobj.(cc).validateInputParams();
            end
            
        end
    end
    
end



%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
