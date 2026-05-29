classdef FittingEIS

    properties
        params0
        Z_re_exp
        Z_im_exp
        omega
        scales
        
        % best_params_found = [] 
        % fitting_error = []

    end



    methods        


        function feis = FittingEIS(params0, scales, Z_re_exp, Z_im_exp, omega)
            
            % feis.params0  = params0;
            % feis.scales   = scales;
            % feis.Z_re_exp = Z_re_exp;
            % feis.Z_im_exp = Z_im_exp;
            % feis.omega    = omega;


            if istable(Z_re_exp), Z_re_exp = table2array(Z_re_exp); end
            if istable(Z_im_exp), Z_im_exp = table2array(Z_im_exp); end
            if istable(omega), omega = table2array(omega); end
            if istable(params0), params0 = table2array(params0); end
            if istable(scales), scales = table2array(scales); end

            % Everything has to be "double" and colon (:)
            feis.params0 = double(params0(:));
            feis.scales  = double(scales(:));
            feis.Z_re_exp= double(Z_re_exp(:));
            feis.Z_im_exp= double(Z_im_exp(:));
            feis.omega   = double(omega(:));

            %pb here because the way I take exp data changes a lot the
            %result...

        end

        function p_norm = unscaled2scaled(feis, p)
            p_norm = (p-feis.params0./feis.scales)./(feis.scales*feis.params0-feis.params0/feis.scales);
        end

        function p = scaled2unscaled(feis, p_norm)
            a = feis.scales;
            p0 = feis.params0;
            p = (a*p0-p0/a).*p_norm +p0/a;
        end


        function [min_value, history, best_params, fitting_error] = optimizationBFGS(feis)

            f_opt = @(p_norm) feis.optifunc(scaled2unscaled(feis, p_norm)) ;
            
            params0_norm = unscaled2scaled(feis, feis.params0);
            % best_params_norm = lsqnonlin(deltagap, params0_norm, lb_norm, ub_norm, [0,0,-1,0,10], 0);
            params0_norm = params0_norm(:);
           

            A_custom = [0, 0, -1, 0, 2*feis.scales(5)/feis.scales(3)];
            b_custom = 0;
            

            [min_value, best_params_norm, history] = unitBoxBFGS(...
                params0_norm, ...                             
                f_opt, ...           
                'maximize', false, ...      
                'linIneq', struct('A', A_custom, 'b', b_custom),  ...   %A*u<=b
                'enforceFeasible', true, ...     
                'maxIt', 300, ...               
                'objChangeTol', 1e-6, ...       
                'gradTol', 1e-5, ...
                'lineSearchMaxIt', 100 ...
                                                                );
            best_params = scaled2unscaled(feis, best_params_norm(:) );

            fitting_error = feis.optifunc(best_params);

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
                g_im_dR1 = -(2*R1*C1.*w) ./ (1+(R1*C1.*w).^2).^2;  
                g_im_dC1 = -(R1^2.*w - C1^2*R1^4.*w.^3) ./ (1+(R1*C1.*w).^2).^2;  
                g_im_dR2 = -(2*R2*C2.*w) ./ (1+(R2*C2.*w).^2).^2;  
                g_im_dC2 = -(R2^2.*w - C2^2*R2^4.*w.^3) ./ (1+(R2*C2.*w).^2).^2; 

                J_im = [g_im_dR0, g_im_dR1, g_im_dC1, g_im_dR2, g_im_dC2];
                derr_re_dp = J_re ./ Modulus_Z; 
                derr_im_dp = J_im ./ Modulus_Z;

                g_true = 2 * (derr_re_dp' * err_re) + 2 * (derr_im_dp' * err_im);
                g_norm = unscaled2scaled(feis, g_true);
                
            end
            

            % finite difference method
            % if nargout > 1
            %     g_norm = zeros(5,1);
            %     dp_norm = 1e-6; 
            % 
            %     for i = 1:5
            %         p_perturb = p;
            %         p_perturb(i) = p_perturb(i) + dp_norm * feis.scales(i);
            % 
            %         [Z_re_pert, Z_im_pert] = load_nyquist(p_perturb, feis.omega);
            % 
            %         err_re_pert = (feis.Z_re_exp(:) - double(Z_re_pert(:))) ./ Modulus_Z;
            %         err_im_pert = (feis.Z_im_exp(:) - double(Z_im_pert(:))) ./ Modulus_Z;
            % 
            %         v_perturb = sum(err_re_pert.^2 + err_im_pert.^2);
            % 
            %         g_norm = (v_perturb - v) / dp_norm;
            %     end
            % end
        end

%% Warburg Model

        function [min_value, history, best_params, fitting_error] = optimizationBFGS_warburg(feis)

            f_opt = @(p_norm) feis.optifunc_warburg(p_norm.* feis.scales) ;
            
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

            fitting_error = feis.optifunc_warburg(best_params);

        end

        function [v, g_norm] = optifunc_warburg(feis, p)
        
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



        %% lsq method

        function residuals = optifunc_lsq(feis, p)
            % 1. On empêche physiquement les paramètres de valoir zéro
            p = max(p, 1e-8); 
            
            % 2. Chargement du modèle
            [Z_re_param, Z_im_param] = load_nyquist(p, feis.omega);
            if istable(Z_re_param), Z_re_param = table2array(Z_re_param); end
            if istable(Z_im_param), Z_im_param = table2array(Z_im_param); end
            
            Z_re_param = double(Z_re_param(:));
            Z_im_param = double(Z_im_param(:));
            
            % Si le modèle renvoie du NaN (physiquement impossible), on pénalise
            if ~isnumeric(Z_re_param) || ~isreal(Z_re_param) || any(isnan(Z_re_param(:)))
                residuals = 1e10 * ones(2 * length(feis.omega), 1);
                return;
            end
            
            % 3. Calcul du module pour l'erreur relative
            Modulus_Z = sqrt(feis.Z_re_exp.^2 + feis.Z_im_exp.^2);
            Modulus_Z = Modulus_Z(:);
            
            % 4. Vecteurs de résidus (Non mis au carré, lsqnonlin s'en charge)
            err_re = (feis.Z_re_exp(:) - Z_re_param(:)) ./ Modulus_Z;
            err_im = (feis.Z_im_exp(:) - Z_im_param(:)) ./ Modulus_Z;
            
            % 5. On empile les parties réelles et imaginaires en un seul vecteur colonne
            residuals = [err_re; err_im];
        end

        function [min_value, best_params, resnorm] = optimizationLsqnonlin(feis)
            % Fonction objectif : renvoie le vecteur de résidus
            f_res = @(p_norm) feis.optifunc_lsq(p_norm(:) .* feis.scales(:));
            
            % Paramètres initiaux normés
            params0_norm = feis.params0(:) ./ feis.scales(:); 
            
            % Bornes (Limites inférieures strictes > 0 pour éviter les crashs)
            lb_norm = ones(5,1) * 1e-5; 
            ub_norm = ones(5,1) * 1e5./ feis.scales(:); 
            
            % Options de lsqnonlin
            options = optimoptions('lsqnonlin', ...
                'Display', 'iter-detailed', ... % Affichage très détaillé
                'Algorithm', 'trust-region-reflective', ... % Parfait pour l'EIS
                'MaxFunctionEvaluations', 3000, ...
                'MaxIterations', 500, ...
                'StepTolerance', 1e-8, ...
                'FunctionTolerance', 1e-6);
                
            % Lancement de l'optimiseur natif MATLAB
            [best_params_norm, resnorm, ~, exitflag, output] = lsqnonlin(...
                f_res, params0_norm, lb_norm, ub_norm, options);
                
            % Récupération des vraies valeurs (dénormalisation)
            best_params = best_params_norm(:) .* feis.scales(:);
            min_value = resnorm; % Erreur finale (somme des carrés des résidus)
            
            % Affichage du statut de fin
            disp('=== RAISON DE FIN DE LSQNONLIN ===');
            disp(output.message);
        end

%% lsq method warburg
        function v = optifunclsq_warburg(feis, p) 
            p_safe = max(p, 1e-15);
            
            R0 = p_safe(1); R1 = p_safe(2); Q1 = p_safe(3); a1 = p_safe(4);
            R2 = p_safe(5); Q2 = p_safe(6); a2 = p_safe(7); Q  = p_safe(8); L = p_safe(9);
            
            w = feis.omega(:);
            s = 1i .* w;
            
            sa1 = s.^a1;
            sa2 = s.^a2;
            Z_param_cplx = R0 + s.*L + R1 ./ (1 + R1.*Q1.*sa1) + R2 ./ (1 + R2.*Q2.*sa2) + 1 ./ (Q.*(s.^0.5));
            
            Z_re_param = real(Z_param_cplx);
            Z_im_param = imag(Z_param_cplx);
            
            if any(isnan(Z_re_param)) || any(isnan(Z_im_param)) || any(isinf(Z_re_param)) || any(isinf(Z_im_param))
                v = 1e6; % Pénalité forte mais gérable
                return;
            end
            
            Modulus_Z = sqrt(Z_re_param.^2 + Z_im_param.^2);
            Modulus_Z = Modulus_Z(:);
            
            err_re = (feis.Z_re_exp(:) - Z_re_param) ./ Modulus_Z;
            err_im = (feis.Z_im_exp(:) - Z_im_param) ./ Modulus_Z;
            
            v = sum(err_re.^2 + err_im.^2);

        end

            
        function [min_value, history, best_params, fitting_error] = optimizationlsq_warburg(feis)
            f_opt = @(p_norm) feis.optifunclsq_warburg(p_norm .* feis.scales);
            
            params0_norm = feis.params0 ./ feis.scales; 
            params0_norm = params0_norm(:);
            
            if length(params0_norm) ~= 9
                error('feis.params0 should have 9 parameters but it has %d).', length(params0_norm));
            end 

            lb_norm = [1e-10; 1e-10; 1e-10; 1e-3; 1e-10; 1e-10; 0.01; 1e-5; 1e-15] ./ feis.scales;
            ub_norm = [1e5;   1e5;   1e6;  1.0;   1e5;  1e6;  1.0;  1e15;   1e-2 ] ./ feis.scales;
            
            options = optimoptions('fmincon', ...
                'Display', 'iter', ...
                'Algorithm', 'sqp', ...       % Très bon pour ce type de fitting
                'SpecifyObjectiveGradient', false, ...
                'MaxFunctionEvaluations', 10000, ... 
                'MaxIterations', 1000); 
                
            [best_params_norm, min_value] = fmincon(f_opt, params0_norm, [], [], [], [], lb_norm, ub_norm, [], options);
            
            history = []; % Si vous n'en avez pas besoin
            best_params = best_params_norm(:) .* feis.scales(:);
            fitting_error = feis.optifunclsq_warburg(best_params);
        end

        %% showing the results

        function plotresults(feis, best_params, fitting_error)

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
            
            % fprintf('\n=== PARAMETERS FOUND ===\n');
            % fprintf('R0 = %.5f Ohms\n', best_params(1));
            % fprintf('R1 = %.5f Ohms\n', best_params(2));
            % fprintf('C1 = %.1f Farads\n', best_params(3));
            % fprintf('R2 = %.5f Ohms\n', best_params(4));
            % fprintf('C2 = %.1f Farads\n', best_params(5));
            
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
