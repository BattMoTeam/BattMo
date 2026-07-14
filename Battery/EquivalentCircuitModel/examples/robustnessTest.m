%% Robustness test for fitting programm
%

for i = 1:5
    
    params0 = pmin + (pmax - pmin) .* rand(1, 5);
    disp(['Tour ', num2str(i), ' - params0 :']);
    disp(params0);
    
    feis = FittingEIS(params0, scales, Z_re_exp, Z_im_exp, omega);

    [~, ~, best_params, fitting_error] = feis.optimizationBFGS();
    feis.plotresults_thevenin(best_params, fitting_error);
    feis.printResults(best_params, fitting_error);

    drawnow; 

end
