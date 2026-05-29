

% filename = 'C:\Users\Alexandre Fichter\Documents\stage_3A\contenu stage\data_August\ank_data\Supplementary material\02_Electrical_characterization\EIS\131-828_EIS_01_MB_CD8.txt';
% filename = '/home/xavier/Matlab/Projects/battmo/Data/131-828_EIS_01_MB_CD8.txt';

% [Z_re_exp, Z_im_exp, omega] = load_experimental_data(filename);
% [Z_re_exp, Z_im_exp, omega] = load_santoni_data();

[Z_re_exp, Z_im_exp, omega] = load_chen_data();


params0 = [3e-2, 1, 3e3, 1e-1, 1e3];
scales = [1e-1, 1e1, 1e5, 1e0, 1e5];

params0_w = [1e-1,1e-1,1e-1,1e-1,1e-2,1,1,1e2,1e-07];
scales_w = params0_w*10;              %to be changed between lsq and bfgs

feis = FittingEIS(params0_w, scales_w, Z_re_exp, Z_im_exp, omega);

[~, ~, best_params, fitting_error] = feis.optimizationBFGS_warburg();

% % Lancement avec LSQNONLIN
% [min_value, best_params, fitting_error] = feis.optimizationLsqnonlin();




%%
feis.plotresults(best_params, fitting_error);

feis.printResults(best_params, fitting_error);



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


% results warburg chen
% 4.7197e-02 ...
% 1.0000e-2 ...
% 2.3392e-03 ...
% 0.0103 ...
% 6.8573e-04 ...
% 1.0000e-1 ...
% 0.0307 ...
% 7.0505e+02 ...
% 2.9134e-06 