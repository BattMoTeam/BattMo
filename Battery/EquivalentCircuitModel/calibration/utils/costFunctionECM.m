function total_error = costFunctionECM(params_test, parameters, initial_values, t_exp, V_exp)


    parameters.R0 = abs(params_test(1));
    parameters.R1 = abs(params_test(2));
    parameters.C1 = max(abs(params_test(3)), 0.1);
    parameters.R2 = abs(params_test(4));
    parameters.C2 = max(abs(params_test(5)), 0.1);

    inputparams = EquivalentCircuitModelInputParams(parameters);
    model = EquivalentCircuitModel(inputparams);
    [t_sim, V_sim, ~] = model.solve();

    % Time synchronization
    V_sim_aligne = interp1(t_sim, V_sim, t_exp, 'linear', 'extrap');

    voltage_error = sum((V_exp - V_sim_aligne).^2);

    lambda = 0.5;

    rel_error = sum((params_test-initial_values)./initial_values);

    tikhonov_error = lambda * rel_error.^2;

    total_error = voltage_error + tikhonov_error;

end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
