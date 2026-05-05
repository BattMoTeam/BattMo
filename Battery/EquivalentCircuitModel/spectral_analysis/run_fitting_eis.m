%% Fitting Script

% initial values [R0, R1, C1, R2, C2]
% params0 = [8.2e-3, 1.5e-3, 1, 1.5e-2, 200];  %  handmade best ones

params0 = [3e-3, 0.00260, 100, 0.00132, 1];


[Z_re_exp, Z_im_exp, omega] = load_experimental_data();


%lengthZ = length(Z_im_exp);

% synthetic data
% opti_params = [0.00341, 0.00260, 5663.9, 0.00132, 0.9];     % params with
% lsqnonlin
% 
% params with unitboxBFGS:
% R0 = 0.00318 Ohms
% R1 = 0.00247 Ohms
% C1 = 10654.4 Farads
% R2 = 0.00155 Ohms
% C2 = 0.4 Farads
% 

% [Z_re_exp, Z_im_exp] = load_nyquist(opti_params, omega);

%%

scales = [10e-2, 10e-3, 15000, 10e-3, 5];

params0_norm = params0 ./ scales; 

lb_true = [1e-5, 1e-5, 1e-5, 1e-5, 1e-5];
lb_norm = lb_true ./ scales;
ub_norm = [];

f_opt = @(p_norm) optifunc(p_norm, scales, Z_re_exp, Z_im_exp, omega) ;

% best_params_norm = lsqnonlin(deltagap, params0_norm, lb_norm, ub_norm, [0,0,-1,0,10], 0);
params0_norm = params0_norm(:);

[min_value, best_params_norm, history] = unitBoxBFGS(...
    params0_norm, ...                             
    f_opt, ...           
    'maximize', false, ...              
    'linIneq', struct('A', [0, 0, -1, 0, 2*scales(5)/scales(3)], 'b', 0), ...      
    'enforceFeasible', true ...         
);


best_params = best_params_norm(:) .* scales(:);

[v_opt, g_opt] = optifunc(best_params, scales, Z_re_exp, Z_im_exp, omega)

fprintf('\n=== FITTING SCORE ===\n');
fprintf('Error : %e\n', v_opt);

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

[Z_re_fit, Z_im_fit] = load_nyquist(best_params, omega);  


%%
figure;
subplot(3,1,1);
semilogx(omega, Z_re_exp, 'ro', 'MarkerFaceColor', 'r');
hold on;
semilogx(omega, Z_re_fit, 'bo', 'LineWidth', 2);        
% legend('Expérimental', 'Modèle (Fitted)');
title('Fitting results');
xlabel('Omega');
ylabel('Z_{re} '); 

subplot(3,1,2);
semilogx(omega, Z_im_exp, 'ro', 'MarkerFaceColor', 'r');
hold on;
semilogx(omega, Z_im_fit, 'bo', 'LineWidth', 2);        
% legend('Expérimental', 'Modèle (Fitted)');
title('Fitting results');
xlabel('Omega');
ylabel('-Z_{im} '); 

subplot(3,1,3);
plot(Z_re_exp, Z_im_exp, 'ro', 'MarkerFaceColor', 'r');
hold on;
plot(Z_re_exp, Z_im_fit, 'bo', 'LineWidth', 2);        
% legend('Expérimental', 'Modèle (Fitted)');
title('Nyquist');
xlabel('Z_{re}');
ylabel('-Z_{im} '); 
axis equal;

grid on;

text_error = sprintf('Fitting error : %.2e', v_opt);

subplot(3,1,3);

text(0.05, 0.90, text_error, 'Units', 'normalized', ...
    'BackgroundColor', 'white', ...   
    'EdgeColor', 'black', ...        
    'FontSize', 11, ...               
    'FontWeight', 'bold');

function [v, g_norm] = optifunc(p_norm, scales, Z_re_exp, Z_im_exp, omega)
     
    p = p_norm(:) .* scales(:);

    [Z_re_param, Z_Im_param] = load_nyquist(p, omega);
    
    Module_Z = sqrt(Z_re_exp.^2 + Z_im_exp.^2);
    Module_Z = Module_Z(:);
    
    err_re = (Z_re_exp(:) - Z_re_param(:)) ./ Module_Z(:);
    err_im = (Z_im_exp(:) - Z_Im_param(:)) ./ Module_Z(:);
    
    v = sum(err_re.^2 + err_im.^2);


if nargout > 1
        w  = omega(:);
        R0 = p(1);
        R1 = p(2);
        C1 = p(3);
        R2 = p(4);
        C2 = p(5);
    
    
        g_re_dR0 = -ones(size(w)); 
        g_re_dR1 = -(1-(R1*C1.*w).^2) ./ (1+(R1*C1.*w).^2).^2;  
        g_re_dC1 = (2*C1*R1^3.*w.^2) ./ (1+(R1*C1.*w).^2).^2;  
        g_re_dR2 = -(1-(R2*C2.*w).^2) ./ (1+(R2*C2.*w).^2).^2;  
        g_re_dC2 = (2*C2*R2^3.*w.^2) ./ (1+(R2*C2.*w).^2).^2;

        J_re = [g_re_dR0, g_re_dR1, g_re_dC1, g_re_dR2, g_re_dC2];
        

        g_im_dR0 = zeros(size(w));
        g_im_dR1 = -(2*R1*C1.*w) ./ (1+(R1*C1.*w).^2).^2;  
        g_im_dC1 = -(R1^2.*w - C1^2*R1^4.*w.^3) ./ (1+(R1*C1.*w).^2).^2;  
        g_im_dR2 = -(2*R2*C2.*w) ./ (1+(R2*C2.*w).^2).^2;  
        g_im_dC2 = -(R2^2.*w - C2^2*R2^4.*w.^3) ./ (1+(R2*C2.*w).^2).^2; 
        
        J_im = [g_im_dR0, g_im_dR1, g_im_dC1, g_im_dR2, g_im_dC2];
        derr_re_dp = J_re ./ Module_Z; 
        derr_im_dp = J_im ./ Module_Z;
    

        g_true = 2 * (derr_re_dp' * err_re) + 2 * (derr_im_dp' * err_im);
        g_norm = g_true .* scales(:);
end

end



