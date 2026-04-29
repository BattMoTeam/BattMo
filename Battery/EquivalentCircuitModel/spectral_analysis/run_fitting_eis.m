%% Fitting Script

% initial values [R0, R1, C1, R2, C2]
 params0 = [8.2e-3, 1.5e-3, 1, 1.5e-2, 200];  %  handmade best ones

% params0 = [6.2261e-3, 0.3e-3, 10000.1264, 0.00097, 2e5];


[Z_re_exp, Z_im_exp, omega] = load_experimental_data();
%lengthZ = length(Z_im_exp);

scales = [1e-3, 1e-3, 1000, 1e-3, 10000];

params0_norm = params0 ./ scales; 

lb_true = [1e-5, 1e-5, 1e-5, 1e-5, 1e-5];
lb_norm = lb_true ./ scales;
ub_norm = [];

delta = @(p_norm) compute_residuals(p_norm .* scales, Z_re_exp, Z_im_exp, omega);

options = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxFunctionEvaluations', 20000);

best_params_norm = lsqnonlin(delta, params0_norm, lb_norm, ub_norm, options);

best_params = best_params_norm .* scales;


fprintf('\n=== PARAMETERS FOUND ===\n');
fprintf('R0 = %.5f Ohms\n', best_params(1));
fprintf('R1 = %.5f Ohms\n', best_params(2));
fprintf('C1 = %.1f Farads\n', best_params(3));
fprintf('R2 = %.5f Ohms\n', best_params(4));
fprintf('C2 = %.1f Farads\n', best_params(5));

%% Graphique 


parameters.R0 = best_params(1);
parameters.R1 = best_params(2);
parameters.C1 = best_params(3);
parameters.R2 = best_params(4);
parameters.C2 = best_params(5);

[Z_re_fit, Z_im_fit] = load_nyquist(best_params, omega);   %  doesn't work with omega which is strange


%%
figure;
subplot(3,1,1);
plot(omega, Z_re_exp, 'ro', 'MarkerFaceColor', 'r');
hold on;
plot(omega, Z_re_fit, 'bo', 'LineWidth', 2);        
legend('Expérimental', 'Modèle (Fitted)');
title('Fitting results');
xlabel('Omega');
ylabel('Z_{re} '); 

subplot(3,1,2);
plot(omega, Z_im_exp, 'ro', 'MarkerFaceColor', 'r');
hold on;
plot(omega, Z_im_fit, 'bo', 'LineWidth', 2);        
legend('Expérimental', 'Modèle (Fitted)');
title('Fitting results');
xlabel('Omega');
ylabel('-Z_{im} '); 

subplot(3,1,3);
plot(Z_re_exp, Z_im_exp, 'ro', 'MarkerFaceColor', 'r');
hold on;
plot(Z_re_exp, Z_im_fit, 'bo', 'LineWidth', 2);        
legend('Expérimental', 'Modèle (Fitted)');
title('Nyquist');
xlabel('Z_{re}');
ylabel('-Z_{im} '); 
axis equal;

grid on;



function res = compute_residuals(p, Z_re_exp, Z_im_exp, omega)

    [Z_re_param, Z_Im_param] = load_nyquist(p, omega);
    
    Module_Z = sqrt(Z_re_exp.^2 + Z_im_exp.^2);
    
    err_re = (Z_re_exp(:) - Z_re_param(:)) ./ Module_Z(:);
    err_im = (Z_im_exp(:) - Z_Im_param(:)) ./ Module_Z(:);
    
    res = [err_re; err_im];
end