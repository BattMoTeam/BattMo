classdef TCtestControlModelInputParams < ControlModelInputParams
%

    properties

        timecontrol

        % for clipping
        lowerCutoffVoltage
        upperCutoffVoltage
        dIdtLimit
        dEdtLimit

        % Tolerances for the control switching (relative tolerances)
        % struct with one value for each switch phase
        % - CC_discharge1 (constant current discharge until lower cutoff voltage is reached)
        % - CC_discharge2 (zero current until dEdtLimit is reached)
        % - CC_charge1    (constant current charge until upper cutoff voltage is reached)
        % - CV_charge2    (constant voltage until dIdtLimit is reached)
        switchTolerances
    end

    methods

        function inputparams = TCtestControlModelInputParams(jsonstruct)

            jsonstruct = setDefaultJsonStructField(jsonstruct, {'switchTolerances', 'CC_discharge1'}, 1e-2);
            jsonstruct = setDefaultJsonStructField(jsonstruct, {'switchTolerances', 'CC_discharge2'}, 0.9);
            jsonstruct = setDefaultJsonStructField(jsonstruct, {'switchTolerances', 'CC_charge1'}, 1e-2);
            jsonstruct = setDefaultJsonStructField(jsonstruct, {'switchTolerances', 'CV_charge2'}, 0.9);

            inputparams = inputparams@ControlModelInputParams(jsonstruct);
            inputparams.controlPolicy = 'TCtest';

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
