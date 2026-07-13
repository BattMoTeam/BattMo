function p = runClassFittingTime(pr)
                                                                
%% Main code:
    time_init = (0:10:1000)';
    current_init = ones(size(time_init));
    
    current_init(time_init >= 0   & time_init < 200)  = 2.5;  % 1. Décharge (0 à 200s)
    current_init(time_init >= 200 & time_init < 400)  = -1.5; % 2. Charge (200 à 400s)
    current_init(time_init >= 400 & time_init < 700)  = 0.0;  % 3. Repos (400 à 700s)
    current_init(time_init >= 700 & time_init <= 1000) = 3.5;  % 4. Décharge (700 à 1000s)

    [time_vec, current_exp, voltage_exp, ocv_vec] = loadP2dVoltage(time_init, current_init);

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
    
    ftime.printResults(best_params, fitting_error);          
    
    ftime.plottest(best_params)
    end

end
