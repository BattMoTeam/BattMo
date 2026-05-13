

% filename = 'C:\Users\Alexandre Fichter\Documents\stage_3A\contenu stage\data_August\ank_data\Supplementary material\02_Electrical_characterization\EIS\131-828_EIS_01_MB_CD8.txt';
% filename = '/home/xavier/Matlab/Projects/battmo/Data/131-828_EIS_01_MB_CD8.txt';
% 
% [Z_re_exp, Z_im_exp, omega] = load_experimental_data(filename);

params0 = [3e-2, 1e-1, 3e1, 1e-2, 1e1];

[Z_re_exp, Z_im_exp, omega] = load_chen_data();

scales = [1e-1, 1, 1e2, 1e-1, 1e2];

feis = FittingEIS(params0, scales, Z_re_exp, Z_im_exp, omega);

% [~, ~, best_params, fitting_error] = feis.optimizationBFGS();

%% Lancement avec LSQNONLIN
[min_value, best_params, fitting_error] = feis.optimizationLsqnonlin();

%%
feis.plotresults(best_params, fitting_error);

feis.printResults(best_params, fitting_error);
