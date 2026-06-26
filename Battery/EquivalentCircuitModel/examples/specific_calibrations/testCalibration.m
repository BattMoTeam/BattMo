omega = logspace(-4, 2, 50);
params = [0.05052, 1.12673, 59119.9, 0.03155, 11054.0];
[Z_re_exp, Z_im_exp] = load_nyquist(params, omega);

params0 = [0.05052, 1.12673, 59119.9, 0.03155, 11054.0];  % initial condition: C1>2*C2

a = 1000;

pmin = params0 / a;
pmax = params0 * a;
scales = [pmin, pmax];

feis = FittingEIS(params0, scales, Z_re_exp, Z_im_exp, omega);

set(0, 'defaultlinelinewidth', 3);
set(0, 'DefaultAxesFontSize', 16);
set(0, 'defaulttextfontsize', 18);

[~, ~, best_params, fitting_error] = feis.optimizationBFGS();
[Z_re_fit, Z_im_fit] = load_nyquist(best_params, feis.omega); 
figure;
plot(feis.Z_re_exp, -feis.Z_im_exp, 'r', 'MarkerFaceColor', 'r');
hold on;
plot(Z_re_fit, -Z_im_fit, 'b');        
legend('experience', 'fitted model');
xlabel('Z_{re}');
ylabel('-Z_{im} '); 
axis equal;
grid on;
text_error = sprintf('Fitting error : %.2e', fitting_error);
text_error = sprintf('Fitting error : %.2e', fitting_error);
text(0.05, 0.85, text_error, 'Units', 'normalized', ...
'BackgroundColor', 'white', ...   
'EdgeColor', 'black', ...        
'FontSize', 11, ...               
'FontWeight', 'bold');
