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

        function inputparams = ElectrodeInputParams(jsonstruct)

            co = 'Coating';
            cc = 'CurrentCollector';


            jsonstruct = equalizeJsonStructField(jsonstruct, {'use_thermal'}, {co, 'use_thermal'});

            include_current_collectors = getJsonStructField(jsonstruct, 'include_current_collectors');
            if isAssigned(include_current_collectors) && include_current_collectors
                jsonstruct = equalizeJsonStructField(jsonstruct, {'use_thermal'}, {cc, 'use_thermal'});
            end

            inputparams = inputparams@ComponentInputParams(jsonstruct);

            pick = @(fd) pickField(jsonstruct, fd);

            inputparams.(co) = CoatingInputParams(pick(co));
            inputparams.(cc) = CurrentCollectorInputParams(pick(cc));

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
