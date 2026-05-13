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
            
            feis.params0  = params0;
            feis.scales   = scales;
            feis.Z_re_exp = Z_re_exp;
            feis.Z_im_exp = Z_im_exp;
            feis.omega    = omega;


            % if istable(Z_re_exp), Z_re_exp = table2array(Z_re_exp); end
            % if istable(Z_im_exp), Z_im_exp = table2array(Z_im_exp); end
            % if istable(omega), omega = table2array(omega); end
            % if istable(params0), params0 = table2array(params0); end
            % if istable(scales), scales = table2array(scales); end
            % 
            % % Everything has to be "double" and colon (:)
            % feis.params0 = double(params0(:));
            % feis.scales  = double(scales(:));
            % feis.Z_re_exp= double(Z_re_exp(:));
            % feis.Z_im_exp= double(Z_im_exp(:));
            % feis.omega   = double(omega(:));

        end


        function [min_value, history, best_params, fitting_error] = optimizationBFGS(feis)

            f_opt = @(p_norm) feis.optifunc(p_norm.* feis.scales) ;
            
            params0_norm = feis.params0 ./ feis.scales; 
            % best_params_norm = lsqnonlin(deltagap, params0_norm, lb_norm, ub_norm, [0,0,-1,0,10], 0);
            params0_norm = params0_norm(:);
           

            A_custom = [0, 0, -1, 0, 2*feis.scales(5)/feis.scales(3)];
            b_custom = 0;
            
            limite_basse = 1e-5;
            A_lim = -eye(5); 
            b_lim = -limite_basse * ones(5, 1);
            
            A_total = [A_custom; A_lim];
            b_total = [b_custom; b_lim];

            [min_value, best_params_norm, history] = unitBoxBFGS(...
                params0_norm, ...                             
                f_opt, ...           
                'maximize', false, ...      
                'linIneq', struct('A', A_total, 'b', b_total), ...   
                'enforceFeasible', true, ...     
                'maxIt', 300, ...               
                'objChangeTol', 1e-6, ...       
                'gradTol', 1e-5, ...
                'lineSearchMaxIt', 100 ...
                                                                );
            best_params = best_params_norm(:) .* feis.scales(:);

            fitting_error = feis.optifunc(best_params);

        end

        function plotresults(feis, best_params, fitting_error)

            [Z_re_fit, Z_im_fit] = load_nyquist(best_params, feis.omega);  
            
            figure;
            subplot(3,1,1);
            semilogx(feis.omega, feis.Z_re_exp, 'r', 'MarkerFaceColor', 'r');
            hold on;
            semilogx(feis.omega, Z_re_fit, 'b', 'LineWidth', 2);        
            legend('experience', 'fitted model');
            title('Fitting results');
            xlabel('Omega');
            ylabel('Z_{re} '); 
            
            subplot(3,1,2);
            semilogx(feis.omega, feis.Z_im_exp, 'r', 'MarkerFaceColor', 'r');
            hold on;
            semilogx(feis.omega, Z_im_fit, 'b', 'LineWidth', 2);        
            legend('experience', 'fitted model');
            title('Fitting results');
            xlabel('Omega');
            ylabel('-Z_{im} '); 
            
            subplot(3,1,3);
            plot(feis.Z_re_exp, feis.Z_im_exp, 'r', 'MarkerFaceColor', 'r');
            hold on;
            plot(feis.Z_re_exp, Z_im_fit, 'b', 'LineWidth', 2);        
            legend('experience', 'fitted model');
            title('Nyquist');
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
            fprintf('R0 = %.5f Ohms\n', best_params(1));
            fprintf('R1 = %.5f Ohms\n', best_params(2));
            fprintf('C1 = %.1f Farads\n', best_params(3));
            fprintf('R2 = %.5f Ohms\n', best_params(4));
            fprintf('C2 = %.1f Farads\n', best_params(5));
            
        end

        function [v, g_norm] = optifunc(feis, p)
            
            [Z_re_param, Z_im_param] = load_nyquist(p, feis.omega);


            if istable(Z_re_param), Z_re_param = table2array(Z_re_param); end
            if istable(Z_im_param), Z_im_param = table2array(Z_im_param); end
            
           
            Z_re_param = double(Z_re_param(:));
            Z_im_param = double(Z_im_param(:));


            if ~isnumeric(Z_re_param) || ~isreal(Z_re_param) || any(isnan(Z_re_param(:)))
                v = 1e10; 
                g_norm = zeros(5,1); 
                return;
            end
            
            Module_Z = sqrt(feis.Z_re_exp.^2 + feis.Z_im_exp.^2);
            Module_Z = Module_Z(:);
            
            err_re = (feis.Z_re_exp(:) - Z_re_param(:)) ./ Module_Z(:);
            err_im = (feis.Z_im_exp(:) - Z_im_param(:)) ./ Module_Z(:);
            
            v = sum(err_re.^2 + err_im.^2);
            
            if nargout > 1

                w  = feis.omega(:);
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
                g_norm = g_true .* feis.scales(:);
                
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
            %         err_re_pert = (feis.Z_re_exp(:) - double(Z_re_pert(:))) ./ Module_Z;
            %         err_im_pert = (feis.Z_im_exp(:) - double(Z_im_pert(:))) ./ Module_Z;
            % 
            %         v_perturb = sum(err_re_pert.^2 + err_im_pert.^2);
            % 
            %         g_norm = (v_perturb - v) / dp_norm;
            %     end
            % end
        end

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
            Module_Z = sqrt(feis.Z_re_exp.^2 + feis.Z_im_exp.^2);
            Module_Z = Module_Z(:);
            
            % 4. Vecteurs de résidus (Non mis au carré, lsqnonlin s'en charge)
            err_re = (feis.Z_re_exp(:) - Z_re_param(:)) ./ Module_Z;
            err_im = (feis.Z_im_exp(:) - Z_im_param(:)) ./ Module_Z;
            
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
            ub_norm = inf(5,1); 
            
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
    end
end
