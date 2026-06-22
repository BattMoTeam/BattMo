filename = 'C:\Users\Alexandre Fichter\Documents\stage_3A\contenu stage\data_August\ank_data\Supplementary material\02_Electrical_characterization\EIS\131-828_EIS_01_MB_CD8.txt';
% filename = '/home/xavier/Matlab/Projects/battmo/Data/131-828_EIS_01_MB_CD8.txt';
[Z_re_exp, Z_im_exp, omega] = load_experimental_data(filename)

params0 = [3.84e-03,...
2.71e-03,...
6.00e+03,...
9.48e-04,...
1.11e+01];              % initial condition: C1>2*C2

a = 1000;

pmin = params0 / a;
pmax = params0 * a;
scales = [pmin, pmax];

feis = FittingEIS(params0, scales, Z_re_exp, Z_im_exp, omega);

[~, ~, best_params, fitting_error] = feis.optimizationBFGS();

[Z_re_fit, Z_im_fit] = load_nyquist(best_params, feis.omega);  
            
figure;

plot(feis.Z_re_exp, -feis.Z_im_exp, 'r', 'MarkerFaceColor', 'r');
hold on;
plot(feis.Z_re_exp, -Z_im_fit, 'b');        
legend('experience', 'fitted model');
xlabel('Z_{re}');
ylabel('-Z_{im} '); 
axis equal;

grid on;

text_error = sprintf('Fitting error : %.2e', fitting_error);

text(0.05, 0.90, text_error, 'Units', 'normalized', ...
     'BackgroundColor', 'white', ...   
     'EdgeColor', 'black', ...        
     'FontSize', 11, ...               
     'FontWeight', 'bold');

fprintf('\n=== FITTING SCORE ===\n');
fprintf('Error : %e\n', fitting_error);

fprintf('\n=== PARAMETERS FOUND ===\n');
fprintf('R0 = %.2e Ohms\n', best_params(1));
fprintf('R1 = %.2e Ohms\n', best_params(2));
fprintf('C1 = %.2e Farads\n', best_params(3));
fprintf('R2 = %.2e Ohms\n', best_params(4));
fprintf('C2 = %.2e Farads\n', best_params(5));