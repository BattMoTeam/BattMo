function p = runClassFittingTime(pr)
                                                                
%% Main code:
    time_init = (0:100:1000)';
    current_init = ones(size(time_init));[time_vec, current_exp, voltage_exp, ocv_vec] = loadP2dVoltage(time_init, current_init);

    params0 = [0.05052, 1.12673, 59119.9, 0.03155, 11054.0];  % initial condition: C1>2*C2
    
    a = 1000;
    
    pmin = params0 / a;
    pmax = params0 * a;
    scales = [pmin, pmax];
        
    ftime = FittingTime(params0, scales, time_vec, current_exp, voltage_exp, ocv_vec);
    
    [~, ~, best_params, fitting_error] = ftime.optimizationBFGS();
    
    p = best_params;
    
    
    %% Printing results
    if nargin == 4 || nargin == 1
    ftime.plotresults_thevenin(best_params, fitting_error);
    
    ftime.printResults(best_params, fitting_error);          % to be changed btw Warburg and Thevenin
    end

end
