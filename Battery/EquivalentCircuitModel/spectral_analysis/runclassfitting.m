params0 = [3e-3, 0.00260, 100, 0.00132, 1];


[Z_re_exp, Z_im_exp, omega] = load_experimental_data();

scales = [10e-2, 10e-3, 15000, 10e-3, 5];

feis = fitting_eis(params0, scales, Z_re_exp, Z_im_exp, omega);

fitting_eis.optimizationBFGS(feis);

fitting_eis.plotresults(feis);

fitting_eis.GiveResults(feis);