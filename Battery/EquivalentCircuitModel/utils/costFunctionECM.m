function total_error = cost_function(params_test, parameters, initial_values, t_exp, V_exp)


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
