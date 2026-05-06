params0 = [3e-3, 0.00260, 100, 0.00132, 1];

filename = 'C:\Users\Alexandre Fichter\Documents\stage_3A\contenu stage\data_August\ank_data\Supplementary material\02_Electrical_characterization\EIS\131-828_EIS_01_MB_CD8.txt';
% filename = '/home/xavier/Matlab/Projects/battmo/Data/131-828_EIS_01_MB_CD8.txt';

[Z_re_exp, Z_im_exp, omega] = load_experimental_data(filename);

scales = [10e-2, 10e-3, 5000, 100e-3, 50];

feis = FittingEIS(params0, scales, Z_re_exp, Z_im_exp, omega);

[~, ~, best_params, fitting_error] = feis.optimizationBFGS();

feis.plotresults(best_params, fitting_error);

feis.printResults(best_params, fitting_error);
