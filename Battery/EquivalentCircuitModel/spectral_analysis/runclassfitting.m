

% filename = 'C:\Users\Alexandre Fichter\Documents\stage_3A\contenu stage\data_August\ank_data\Supplementary material\02_Electrical_characterization\EIS\131-828_EIS_01_MB_CD8.txt';
% filename = '/home/xavier/Matlab/Projects/battmo/Data/131-828_EIS_01_MB_CD8.txt';
% 
% [Z_re_exp, Z_im_exp, omega] = load_experimental_data(filename);

params0 = [3e-2, 1, 3e3, 1e-1, 1e3];

[Z_re_exp, Z_im_exp, omega] = load_chen_data();

scales = [1e-1, 1e1, 1e5, 1e0, 1e5];

% feis = FittingEIS(params0, scales, Z_re_exp, Z_im_exp, omega);
% 
% [~, ~, best_params, fitting_error] = feis.optimizationBFGS();

% Lancement avec LSQNONLIN
[min_value, best_params, fitting_error] = feis.optimizationLsqnonlin();


% result with lsqnonlin:
% === FITTING SCORE ===
% Error : 7.324871e-02
% 
% === PARAMETERS FOUND ===
% R0 = 0.05052 Ohms
% R1 = 1.12673 Ohms
% C1 = 59119.9 Farads
% R2 = 0.03155 Ohms
% C2 = 11054.0 Farads


%%
feis.plotresults(best_params, fitting_error);

feis.printResults(best_params, fitting_error);


