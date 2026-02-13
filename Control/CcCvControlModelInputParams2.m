classdef CcCvControlModelInputParams2 < ControlModelInputParams
%
% Constant-Current-Constant-Voltage Control. We switch between charge and discharge scenarii. The switch occurs when the
% derivative of the non-controlled variable gets below a given threshold
%

    properties

        times
        values

        %
        % Charge and discharge rates
        CRate
        DRate

        lowerCutoffVoltage
        upperCutoffVoltage

        % When voltage control, we wait for the derivative of the current to be below the value of dIdtLimit between we
        % switch to the following constant current control.
        dIdtLimit

        % When we switch to current control equal to zero after discharge, we wait for the voltage to stabilize and
        % switch to charge when the derivative of the voltage with respect to time is lower than the threshold given by
        % dEdtLimit.
        dEdtLimit

        % Control used initially. String that can take one of the following values
        % - 'discharging'
        % - 'charging'
        initialControl

        % Number of cycles
        numberOfCycles

        % Tolerances for the control switching (relative tolerances)
        % struct with one value for each switch phase
        % - CC_discharge1 (constant current discharge until lower cutoff voltage is reached)
        % - CC_discharge2 (zero current until dEdtLimit is reached)
        % - CC_charge1    (constant current charge until upper cutoff voltage is reached)
        % - CV_charge2    (constant voltage until dIdtLimit is reached)
        switchTolerances
    end

    methods

        function inputparams = CcCvControlModelInputParams2(jsonstruct)


            jsonstruct = setDefaultJsonStructField(jsonstruct, {'switchTolerances', 'CC_discharge1'}, 1e-2);
            jsonstruct = setDefaultJsonStructField(jsonstruct, {'switchTolerances', 'CC_discharge2'}, 0.9);
            jsonstruct = setDefaultJsonStructField(jsonstruct, {'switchTolerances', 'CC_charge1'}, 1e-2);
            jsonstruct = setDefaultJsonStructField(jsonstruct, {'switchTolerances', 'CV_charge2'}, 0.9);

            inputparams = inputparams@ControlModelInputParams(jsonstruct);
            inputparams.controlPolicy = 'CCCV2';

            inputparams.times = jsonstruct.times;
            inputparams.values = jsonstruct.values;

            assert(all(isfinite(inputparams.times)), 'Times must be finite');
            assert(all(diff(inputparams.times) > 0), 'Times must be strictly increasing');
            assert(all(isfinite(inputparams.values)), 'Values must be finite');

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
