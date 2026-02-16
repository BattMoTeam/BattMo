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

        %% coating model setup structure

        coatingModelSetup % structure that determines the type of coating with field
                          % - swelling : boolean, true if swelling coating is used

        %% Parameters assigned at setup

        include_current_collectors
        use_thermal
        use_normed_current_collector
        
    end

    methods

        function inputparams = ElectrodeInputParams(jsonstruct)

            co = 'Coating';
            am = 'ActiveMaterial';
            cc = 'CurrentCollector';

            jsonstruct = equalizeJsonStructField(jsonstruct, {'use_thermal'}, {co, 'use_thermal'});

            include_current_collectors = getJsonStructField(jsonstruct, 'include_current_collectors');

            if isAssigned(include_current_collectors) && include_current_collectors
                jsonstruct = equalizeJsonStructField(jsonstruct, {'use_thermal'}, {cc, 'use_thermal'});
            end

            if include_current_collectors
                jsonstruct = setDefaultJsonStructField(jsonstruct, {'use_normed_current_collector'}, false);
            end

            jsonstruct = setDefaultJsonStructField(jsonstruct, {'coatingModelSetup', 'swelling'}, false);
            
            jsonstruct = equalizeJsonStructField(jsonstruct, ...
                                                 {'coatingModelSetup', 'swelling'}, ...
                                                 {'Coating', 'activeMaterialModelSetup', 'swelling'});
            
            is_swelling = getJsonStructField(jsonstruct, {'coatingModelSetup', 'swelling'});

            if is_swelling
                jsonstruct = setJsonStructField(jsonstruct, {co, am, 'diffusionModelType'}, 'swelling');
            end
            
            inputparams = inputparams@ComponentInputParams(jsonstruct);

            pick = @(fd) pickField(jsonstruct, fd);

            if is_swelling
                inputparams.(co) = SwellingCoatingInputParams(pick(co));
            else
                inputparams.(co) = CoatingInputParams(pick(co));
            end
            
            if include_current_collectors
                if getJsonStructField(jsonstruct, 'use_normed_current_collector')
                    inputparams.(cc) = NormedCurrentCollectorInputParams(pick(cc));
                else
                    inputparams.(cc) = CurrentCollectorInputParams(pick(cc));
                end
            end

        end

    end

end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
