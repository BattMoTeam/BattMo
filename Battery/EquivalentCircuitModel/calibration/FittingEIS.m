classdef FittingEIS < handle

    properties

        % Observed impedances (array of same lengths)
        omega 
        Z_re_exp
        Z_im_exp

        % initial parameters
        params0

        % scales for unit box search algorithm
        scales
        
    end



    methods        


        function feis = FittingEIS(eisdata, params0, scales)
            
            feis.Z_re_exp = eisdata.Z_re_exp;
            feis.Z_im_exp = eisdata.Z_im_exp;
            feis.omega    = eisdata.omega;
            feis.params0  = params0;
            
            if nargin < 3 || isempty(scales)
                feis.setupDefaultScales();
            end

        end

        function setupDefaultScales(feis)

            params0 = feis.params0;

            a = 1000;
            
            pmin = params0 / a;
            pmax = params0 * a;
            scales = [pmin; pmax];

            feis.scales = scales;
            
        end
        
        function p_norm = unscaled2scaled(feis, p)
            
            pmin = feis.scales(1:5);
            pmax = feis.scales(6:10);
            p_norm = (p-pmin)./(pmax-pmin);
            
        end

        function p = scaled2unscaled(feis, p_norm)
            
            pmin = feis.scales(1:5);
            pmax = feis.scales(6:10);
            p = (pmax-pmin).*p_norm +pmin;
            
        end


        function [optparams, fitting_error] = run(feis)
            
            [~, ~, optparams, fitting_error] = feis.optimizationBFGS();
            
        end
        
        %% Thevenin Model

        function [min_value, history, best_params, fitting_error] = optimizationBFGS(feis)

            f_opt = @(p_norm) feis.optifunc(scaled2unscaled(feis, p_norm));
            
            params0_norm = feis.unscaled2scaled(feis.params0);
            % best_params_norm = lsqnonlin(deltagap, params0_norm, lb_norm, ub_norm, [0,0,-1,0,10], 0);
            params0_norm = params0_norm(:);
           
            pmin = feis.scales(1:5);
            pmax = feis.scales(6:10);
            p03min = pmin(3);
            p03max = pmax(3);
            p05min = pmin(5);
            p05max = pmax(5);

            A_custom = [0, 0, -1, 0, 2*(p05max - p05min)/(p03max- p03min)];
            b_custom = (p03min-2*p05min)/(p03max - p03min);
            
            max_iter = 300;
            tol_obj = 1e-15;
            tol_grad = 1e-5;

            [min_value, best_params_norm, history] = unitBoxBFGS(...
                params0_norm, ...                             
                f_opt, ...           
                'maximize'       , false, ...      
                'linIneq'        , struct('A', A_custom, 'b', b_custom), ... 
                'enforceFeasible', true    , ...     
                'maxIt'          , max_iter, ...               
                'objChangeTol'   , tol_obj , ...       
                'gradTol'        , tol_grad, ...
                'lineSearchMaxIt', 100);
            best_params = scaled2unscaled(feis, best_params_norm(:) );

            % explications of why it stopped
            fitting_error = feis.optifunc(best_params);

            it_count = length(history.val) - 1; 
            final_pg = history.pg(end);         
            
            % Cost function variation
            if it_count >= 1
                delta_v = abs(history.val(end) - history.val(end-1));
            else
                delta_v = inf;
            end
            
            if it_count >= max_iter
                fprintf('Stopped because %d iterations reached.\n', max_iter);
            elseif final_pg < tol_grad
                fprintf('Stopped because reached gradTol = %e.\n', tol_grad);
            elseif delta_v < tol_obj
                fprintf('Stopped because reached objChangeTol = %e.\n', tol_obj);
            else
                fprintf('Unexpected error or relative tolerance reached.\n');
            end           

        end

        function [v, g_norm] = optifunc(feis, p)
            
            % Preventing physical values to get too close from 0
            p_safe = max(p, 1e-8);

            [Z_re_param, Z_im_param] = load_nyquist(p_safe, feis.omega);


            if istable(Z_re_param), Z_re_param = table2array(Z_re_param); end
            if istable(Z_im_param), Z_im_param = table2array(Z_im_param); end
            
           
            Z_re_param = double(Z_re_param(:));
            Z_im_param = double(Z_im_param(:));


            if ~isnumeric(Z_re_param) || ~isreal(Z_re_param) || any(isnan(Z_re_param(:)))
                v = 1e10; 
                g_norm = zeros(5,1); 
                return;
            end
            
            Modulus_Z = sqrt(feis.Z_re_exp.^2 + feis.Z_im_exp.^2);
            Modulus_Z = Modulus_Z(:);
            
            err_re = (feis.Z_re_exp(:) - Z_re_param(:)) ./ Modulus_Z(:);
            err_im = (feis.Z_im_exp(:) - Z_im_param(:)) ./ Modulus_Z(:);
            
            v = sum(err_re.^2 + err_im.^2);


            if nargout > 1

                w  = feis.omega(:);
                R0 = p_safe(1);
                R1 = p_safe(2);
                C1 = p_safe(3);
                R2 = p_safe(4);
                C2 = p_safe(5);

                g_re_dR0 = -ones(size(w)); 
                g_re_dR1 = -(1-(R1*C1.*w).^2) ./ (1+(R1*C1.*w).^2).^2;  
                g_re_dC1 = (2*C1*R1^3.*w.^2) ./ (1+(R1*C1.*w).^2).^2;  
                g_re_dR2 = -(1-(R2*C2.*w).^2) ./ (1+(R2*C2.*w).^2).^2;  
                g_re_dC2 = (2*C2*R2^3.*w.^2) ./ (1+(R2*C2.*w).^2).^2;

                J_re = [g_re_dR0, g_re_dR1, g_re_dC1, g_re_dR2, g_re_dC2];

                g_im_dR0 = zeros(size(w));
                g_im_dR1 = (2*R1*C1.*w) ./ (1+(R1*C1.*w).^2).^2;  
                g_im_dC1 = (R1^2.*w - C1^2*R1^4.*w.^3) ./ (1+(R1*C1.*w).^2).^2;  
                g_im_dR2 = (2*R2*C2.*w) ./ (1+(R2*C2.*w).^2).^2;  
                g_im_dC2 = (R2^2.*w - C2^2*R2^4.*w.^3) ./ (1+(R2*C2.*w).^2).^2; 

                J_im = [g_im_dR0, g_im_dR1, g_im_dC1, g_im_dR2, g_im_dC2];
                derr_re_dp = J_re ./ Modulus_Z; 
                derr_im_dp = J_im ./ Modulus_Z;

                g_true = 2 * (derr_re_dp' * err_re) + 2 * (derr_im_dp' * err_im);
                % Chain rule : 
                pmin = feis.scales(1:5);
                pmax = feis.scales(6:10);
                dp_dpnorm = (pmax - pmin);      % Dérivée de p par rapport à p_norm
                g_norm = g_true .* dp_dpnorm;
               
                
            end
        end 

        %% showing the results

        function plotResults(feis, best_params, fitting_error)

            [Z_re_fit, Z_im_fit] = load_nyquist(best_params, feis.omega);  
            
            figure;
            subplot(3,1,1);
            semilogx(feis.omega, feis.Z_re_exp, 'r', 'MarkerFaceColor', 'r');
            hold on;
            semilogx(feis.omega, Z_re_fit, 'b');        
            legend('experience', 'fitted model');
            title('Fitting Real Impedance');
            xlabel('Omega');
            ylabel('Z_{re} '); 
            
            subplot(3,1,2);
            semilogx(feis.omega, -feis.Z_im_exp, 'r', 'MarkerFaceColor', 'r');
            hold on;
            semilogx(feis.omega, -Z_im_fit, 'b');        
            legend('experience', 'fitted model');
            title('Fitting Imaginary Impedance');
            xlabel('Omega');
            ylabel('-Z_{im} '); 
            
            subplot(3,1,3);
            plot(feis.Z_re_exp, -feis.Z_im_exp, 'r', 'MarkerFaceColor', 'r');
            hold on;
            plot(feis.Z_re_exp, -Z_im_fit, 'b');        
            legend('experience', 'fitted model');
            title('Nyquist Diagram');
            xlabel('Z_{re}');
            ylabel('-Z_{im} '); 
            axis equal;
            
            grid on;
            
            text_error = sprintf('Fitting error : %.2e', fitting_error);
            
            subplot(3,1,3);
            
            text(0.05, 0.90, text_error, 'Units', 'normalized', ...
                 'BackgroundColor', 'white', ...   
                 'EdgeColor', 'black', ...        
                 'FontSize', 11, ...               
                 'FontWeight', 'bold');
        end


        function printResults(feis, best_params, fitting_error)

            fprintf('\n=== FITTING SCORE ===\n');
            fprintf('Error : %e\n', fitting_error);
            
            fprintf('\n=== PARAMETERS FOUND ===\n');
            fprintf('R0 = %.2e Ohms\n', best_params(1));
            fprintf('R1 = %.2e Ohms\n', best_params(2));
            fprintf('C1 = %.2e Farads\n', best_params(3));
            fprintf('R2 = %.2e Ohms\n', best_params(4));
            fprintf('C2 = %.2e Farads\n', best_params(5));
            
        end
        
    end
end
