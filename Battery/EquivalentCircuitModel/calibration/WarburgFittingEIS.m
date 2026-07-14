classdef  WarburgFittingEIS < FittingEIS

    methods        

        function feis = WarburgFittingEIS(params0, eisdata, scales)

            if nargin < 3
                scales = [];
            end
            
            feis = feis@FittingEIS(params0, eisdata, scales);

        end

        function [min_value, history, best_params, fitting_error] = optimizationBFGS(feis)

            %% Warburg Model
            error('not updated....');
            
            f_opt = @(p_norm) feis.optifunc(p_norm.* feis.scales) ;
            
            params0_norm = feis.params0 ./ feis.scales; 
            params0_norm = params0_norm(:);
          
            

            [min_value, best_params_norm, history] = unitBoxBFGS(...
                params0_norm, ...                             
                f_opt, ...           
                'maximize', false, ...      
                'enforceFeasible', true, ...     
                'maxIt', 300, ...               
                'objChangeTol', 1e-6, ...       
                'gradTol', 1e-5, ...
                'lineSearchMaxIt', 200 ...
                                           );
            best_params = best_params_norm(:) .* feis.scales(:);

            fitting_error = feis.optifunc(best_params);

        end

        function [v, g_norm] = optifunc(feis, p)
        
            % Preventing physical values to get too close from 0
            p_safe = max(p, 1e-7);


            w  = feis.omega(:);
            R0 = p_safe(1);
            R1 = p_safe(2);
            Q1 = p_safe(3);
            a1 = p_safe(4);
            R2 = p_safe(5);
            Q2 = p_safe(6);
            a2 = p_safe(7);
            Q  = p_safe(8);
            L  = p_safe(9);

            s = 1i .* w;
            s = s(:);
            
            sa1 = s.^a1;
            sa2 = s.^a2;

            Z_param_cplx = R0 + s.*L + R1 ./ (1 + R1.*Q1.*sa1) + R2 ./ (1 + R2.*Q2.*sa2) + 1 ./ (Q.*(s.^0.5));
    
            Z_re_param = real(Z_param_cplx);
            Z_im_param = imag(Z_param_cplx);

            if istable(Z_re_param), Z_re_param = table2array(Z_re_param); end
            if istable(Z_im_param), Z_im_param = table2array(Z_im_param); end
            
           
            Z_re_param = double(Z_re_param(:));
            Z_im_param = double(Z_im_param(:));

            Modulus_Z = sqrt(Z_re_param.^2 + Z_im_param.^2);
            Modulus_Z = Modulus_Z(:);
            
            err_re = (feis.Z_re_exp(:) - Z_re_param(:)) ./ Modulus_Z(:);
            err_im = (feis.Z_im_exp(:) - Z_im_param(:)) ./ Modulus_Z(:);
            
            v = sum(err_re.^2 + err_im.^2);
             
            
            if nargout > 1

            den1 = (1 + R1 .* Q1 .* sa1).^2;
            den2 = (1 + R2 .* Q2 .* sa2).^2;
            
            
            dZ_dR0 = ones(size(s));                     %  R0
            dZ_dR1 = 1 ./ den1;                         %  R1
            dZ_dR2 = 1 ./ den2;                         %  R2
            
            dZ_dQ1 = - (R1^2 .* sa1) ./ den1;           % Q1
            dZ_dQ2 = - (R2^2 .* sa2) ./ den2;           % Q2
            
            dZ_da1 = - (Q1 * R1^2 .* sa1 .* log(s)) ./ den1; % a1
            dZ_da2 = - (Q2 * R2^2 .* sa2 .* log(s)) ./ den2; % a2
            
            dZ_dQ  = - 1 ./ (Q^2 .* (s.^0.5));          % Q
            dZ_dL  = s;                                 % L
            
            J = [dZ_dR0, dZ_dR1, dZ_dQ1, dZ_da1, dZ_dR2, dZ_dQ2, dZ_da2, dZ_dQ, dZ_dL];

            J_im = imag(J);
            J_re = real(J);

            derr_re_dp = -J_re ./ Modulus_Z; 
            derr_im_dp = -J_im ./ Modulus_Z;

            g_true = 2 * (derr_re_dp' * err_re) + 2 * (derr_im_dp' * err_im);
            g_norm = g_true .* feis.scales(:);
            
            
            if any(isinf(g_norm)) || any(isnan(g_norm)) || isinf(v) || isnan(v)  v = 1e10; 
                g_norm = zeros(9,1); 
                return;
            end
            end
        end

        function plotResults(feis, best_params, fitting_error)

            % [Z_re_fit, Z_im_fit] = load_nyquist(best_params, feis.omega);  
            w  = feis.omega(:);
            R0 = best_params(1);
            R1 = best_params(2);
            Q1 = best_params(3);
            a1 = best_params(4);
            R2 = best_params(5);
            Q2 = best_params(6);
            a2 = best_params(7);
            Q  = best_params(8);
            L  = best_params(9);
            s = 1i .* w;
            s = s(:);
            
            sa1 = s.^a1;
            sa2 = s.^a2;

            Z_param_cplx = R0 + s.*L + R1 ./ (1 + R1.*Q1.*sa1) + R2 ./ (1 + R2.*Q2.*sa2) + 1 ./ (Q.*(s.^0.5));
    
            Z_re_fit = real(Z_param_cplx);
            Z_im_fit = imag(Z_param_cplx);


            figure;
            subplot(3,1,1);
            semilogx(feis.omega, feis.Z_re_exp, 'r', 'MarkerFaceColor', 'r');
            hold on;
            semilogx(feis.omega, Z_re_fit, 'b');        
            legend('experience', 'fitted model');
            title('Fitting results');
            xlabel('Omega');
            ylabel('Z_{re} '); 
            
            subplot(3,1,2);
            semilogx(feis.omega, -feis.Z_im_exp, 'r', 'MarkerFaceColor', 'r');
            hold on;
            semilogx(feis.omega, -Z_im_fit, 'b');        
            legend('experience', 'fitted model');
            title('Fitting results');
            xlabel('Omega');
            ylabel('-Z_{im} '); 
            
            subplot(3,1,3);
            plot(feis.Z_re_exp, -feis.Z_im_exp, 'r', 'MarkerFaceColor', 'r');
            hold on;
            plot(Z_re_fit, -Z_im_fit, 'b');        
            legend('experience', 'fitted model');
            title('Nyquist');
            xlabel('Z_{re}');
            ylabel('-Z_{im} '); 
            axis equal;
            
            grid on;
            text_error = sprintf('Fitting error : %.2e', fitting_error);
            text(0.05, 0.85, text_error, 'Units', 'normalized', ...
            'BackgroundColor', 'white', ...   
            'EdgeColor', 'black', ...        
            'FontSize', 11, ...               
            'FontWeight', 'bold');
        end

        function printResults(feis, best_params, fitting_error)

            fprintf('\n=== FITTING SCORE ===\n');
            fprintf('Error : %e\n', fitting_error);
            

            fprintf('\n=== PARAMETERS FOUND ===\n');
            fprintf('R0    = %.4e Ohms\n', best_params(1));
            fprintf('R1    = %.4e Ohms\n', best_params(2));
            fprintf('Q1    = %.4e s^a/Ohm\n', best_params(3));
            fprintf('a1    = %.4f \n', best_params(4)); % %f suffit pour 'a' car il est entre 0 et 1
            fprintf('R2    = %.4e Ohms\n', best_params(5));
            fprintf('Q2    = %.4e s^a/Ohm\n', best_params(6));
            fprintf('a2    = %.4f \n', best_params(7));
            fprintf('Q     = %.4e \n', best_params(8));
            fprintf('L     = %.4e Henrys\n', best_params(9));
            
        end
    end
end
