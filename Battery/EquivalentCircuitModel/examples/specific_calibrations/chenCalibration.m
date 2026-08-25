%% ECM calibration with P2D model
%
% We calibrate the parameters of the equivalent circuit model (ECM) using the P2D model. The
% calibration is performed by fitting the impedance data obtained from the P2D model to the ECM.
%

%%
% clear workspace

clear all

%%
% In the function :battmo:`load_chen_data`, we load the impedance data from the P2D model. The data
% is stored in a structure called `eisdata`, which contains the frequencies and the real and
% imaginary parts of the Impedances. The P2D model is set using the Chen et al (2020) dataset 

eisdata = load_chen_data();

%%
% Intial guess for the ECM parameters. The value corresponds to the five parameters of the ECM
% model, given in the following order, $R0$, $R1$, $C1$, $R2$ and $C2$.
%
params0 = [5e-2; 1; 6e4; 3e-2; 1e1];  

%%
% We use a scaling factor to define lower and upper parameter bounds (pmin and pmax) that are used
% in the optimization algorithm. These bounds restrict the search domain and can be used to speed up
% convergence.
%
a = 1e1;
pmin = params0 / a;
pmax = params0 * a;
scales = [pmin; pmax];

%% Optimization
% We setup and run the EIS fitting algorithm
%
feis = FittingEIS(eisdata, params0, scales);

[~, ~, best_params, fitting_error] = feis.optimizationBFGS();

%%
%  We compute the impedance values form the ECM parameters using the function :battmo:`load_nyquist`
%
[Z_re_fit, Z_im_fit] = load_nyquist(best_params, feis.omega);  

%% Plotting
% 

set(0, 'defaultlinelinewidth', 3);
set(0, 'DefaultAxesFontSize', 16);
set(0, 'defaulttextfontsize', 18);


figure;

plot(feis.Z_re_exp, -feis.Z_im_exp, 'r', 'MarkerFaceColor', 'r');
hold on;
plot(feis.Z_re_exp, -Z_im_fit, 'b');        
legend('experience', 'fitted model');
xlabel('Z_{re}');
ylabel('-Z_{im} '); 
axis equal;

grid on;

