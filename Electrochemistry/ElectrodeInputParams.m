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

            inputparams = inputparams@ComponentInputParams(jsonstruct);

            co = 'Coating';
            cc = 'CurrentCollector';

            pick = @(fd) pickField(jsonstruct, fd);

            inputparams.(co) = CoatingInputParams(pick(co));
            inputparams.(cc) = CurrentCollectorInputParams(pick(cc));

            inputparams = inputparams.validateInputParams();

        end

        function inputparams = validateInputParams(inputparams)

            am = 'Coating';
            cc = 'CurrentCollector';

            inputparams = mergeParameters(inputparams, {{'use_thermal'}, {am, 'use_thermal'}});
            inputparams.(am) = inputparams.(am).validateInputParams();

            if ~isempty(inputparams.include_current_collectors)
                inputparams = mergeParameters(inputparams, {{'use_thermal'}, {cc, 'use_thermal'}});
                inputparams.(cc) = inputparams.(cc).validateInputParams();
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
