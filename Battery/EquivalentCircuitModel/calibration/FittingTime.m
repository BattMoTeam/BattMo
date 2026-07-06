classdef FittingTime

    properties
        params0
        scales
        time_vec
        current_exp
        voltage_exp
        ocv_vec
    end


    methods        
    
        function ftime = FittingTime(params0, scales, time_vec, current_exp, voltage_exp, ocv_vec)
           

            if istable(params0), params0 = table2array(params0); end
            if istable(scales), scales = table2array(scales); end
            if istable(time_vec), time_vec = table2array(time_vec); end
            if istable(current_exp), current_exp = table2array(current_exp); end
            if istable(voltage_exp), voltage_exp = table2array(voltage_exp); end
            if istable(ocv_vec), ocv_vec = table2array(ocv_vec); end
            
            ftime.params0     = double(params0(:));
            ftime.scales      = double(scales(:));
            ftime.time_vec    = double(time_vec(:));
            ftime.current_exp = double(current_exp(:));
            ftime.voltage_exp = double(voltage_exp(:));
            ftime.ocv_vec     = double(ocv_vec(:));
        end


        function p_norm = unscaled2scaled(ftime, p)
            pmin = ftime.scales(1:5);
            pmax = ftime.scales(6:10);
            p_norm = (p-pmin)./(pmax-pmin);
        end

        function p = scaled2unscaled(ftime, p_norm)
            pmin = ftime.scales(1:5);
            pmax = ftime.scales(6:10);
            p = (pmax-pmin).*p_norm +pmin;
        end

        
        function [min_value, history, best_params, fitting_error] = optimizationBFGS(ftime)
            

            f_opt = @(p_norm) ftime.optifunc(scaled2unscaled(ftime, p_norm));
            
            params0_norm = unscaled2scaled(ftime, ftime.params0);
            % best_params_norm = lsqnonlin(deltagap, params0_norm, lb_norm, ub_norm, [0,0,-1,0,10], 0);
            params0_norm = params0_norm(:);
           
            pmin = ftime.scales(1:5);
            pmax = ftime.scales(6:10);
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
                'maximize', false, ...      
                'linIneq', struct('A', A_custom, 'b', b_custom),  ...   %A*u<=b
                'enforceFeasible', true, ...     
                'maxIt', max_iter, ...               
                'objChangeTol', tol_obj, ...       
                'gradTol', tol_grad, ...
                'lineSearchMaxIt', 100 ...
                                                                );
            best_params = scaled2unscaled(ftime, best_params_norm(:) );

            % explications of why it stopped
            fitting_error = ftime.optifunc(best_params);

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

        function [v, g_norm] = optifunc(ftime, p)
            
            % Preventing physical values to get too close from 0
            p_safe = max(p, 1e-8);
            R0 = p_safe(1);
            R1 = p_safe(2);
            C1 = p_safe(3);
            R2 = p_safe(4);
            C2 = p_safe(5);
            
            N_points = length(ftime.time_vec);
            V_ecm = zeros(N_points, 1);
            

            Uc1 = 0; 
            Uc2 = 0;
            
            % Sensibilities initialisation 
            
            s_R1 = 0; % dUc1 / dR1
            s_C1 = 0; % dUc1 / dC1
            s_R2 = 0; % dUc2 / dR2
            s_C2 = 0; % dUc2 / dC2
            
            g_true = zeros(5, 1);
            
            for k = 1:N_points
                if k == 1
                    dt = ftime.time_vec(1);
                else
                    dt = ftime.time_vec(k) - ftime.time_vec(k-1);
                end
                
                I = ftime.current_exp(k);           
                V_ocv = ftime.ocv_vec(k); 
                
                % 1. Tension du modèle au pas k
                V_ecm(k) = V_ocv - R0 * I - Uc1 - Uc2;
                
                % 2. Erreur instantanée
                err = ftime.voltage_exp(k) - V_ecm(k);
                
                % 3. Si le gradient est demandé (nargout > 1)
                if nargout > 1
                    % Dérivées partielles de V_ecm par rapport à chaque paramètre
                    dV_dR0 = -I;
                    dV_dR1 = -s_R1;
                    dV_dC1 = -s_C1;
                    dV_dR2 = -s_R2;
                    dV_dC2 = -s_C2;
                    
                    % Accumulation du gradient de la fonction coût (v = sum(err^2))
                    g_true(1) = g_true(1) - 2 * err * dV_dR0;
                    g_true(2) = g_true(2) - 2 * err * dV_dR1;
                    g_true(3) = g_true(3) - 2 * err * dV_dC1;
                    g_true(4) = g_true(4) - 2 * err * dV_dR2;
                    g_true(5) = g_true(5) - 2 * err * dV_dC2;
                    
                    % Mise à jour des sensibilités pour le pas suivant (Euler Explicite)
                    s_R1 = s_R1 + (-s_R1 / (R1 * C1) + Uc1 / (R1^2 * C1)) * dt;
                    s_C1 = s_C1 + (-s_C1 / (R1 * C1) + Uc1 / (R1 * C1^2) - I / C1^2) * dt;
                    s_R2 = s_R2 + (-s_R2 / (R2 * C2) + Uc2 / (R2^2 * C2)) * dt;
                    s_C2 = s_C2 + (-s_C2 / (R2 * C2) + Uc2 / (R2 * C2^2) - I / C2^2) * dt;
                end
                
                % 4. Mise à jour des variables d'état électriques pour le pas suivant
                Uc1 = Uc1 + (-Uc1 / (R1 * C1) + I / C1) * dt;
                Uc2 = Uc2 + (-Uc2 / (R2 * C2) + I / C2) * dt;
            end
            
            % Valeur de la fonction objectif (Erreur quadratique totale)
            v = sum((ftime.voltage_exp - V_ecm).^2);
            
            % 5. Application de la règle de dérivation en chaîne pour l'espace normalisé
            if nargout > 1
                pmin = ftime.scales(1:5);
                pmax = ftime.scales(6:10);
                dp_dpnorm = (pmax - pmin);      
                g_norm = g_true .* dp_dpnorm(:);
            end
        end 
    
        function plotresults_thevenin(ftime, best_params, fitting_error)
            % 1. Recalculate the final voltage trajectory using best_params
            p_safe = max(best_params, 1e-8);
            R0 = p_safe(1);
            R1 = p_safe(2);
            C1 = p_safe(3);
            R2 = p_safe(4);
            C2 = p_safe(5);
            
            N_points = length(ftime.time_vec);
            V_ecm = zeros(N_points, 1);
            Uc1 = 0; 
            Uc2 = 0;
            
            for k = 1:N_points
                if k == 1
                    dt = ftime.time_vec(1);
                else
                    dt = ftime.time_vec(k) - ftime.time_vec(k-1);
                end
                I = ftime.current_exp(k);           
                V_ocv = ftime.ocv_vec(k); 
                
                V_ecm(k) = V_ocv - R0 * I - Uc1 - Uc2;
                
                Uc1 = Uc1 + (-Uc1 / (R1 * C1) + I / C1) * dt;
                Uc2 = Uc2 + (-Uc2 / (R2 * C2) + I / C2) * dt;
            end
            
            % 2. Generate time-domain plots
            figure('Name', 'ECM Time Domain Fitting Results', 'NumberTitle', 'off');
            
            % Subplot 1: Applied Current Profile
            subplot(3,1,1);
            plot(ftime.time_vec, ftime.current_exp, 'k', 'LineWidth', 1.5);
            grid on;
            title('Applied Current Profile');
            xlabel('Time (s)');
            ylabel('Current (A)');
            
            % Subplot 2: Voltage Comparison (Exp vs Model)
            subplot(3,1,2);
            plot(ftime.time_vec, ftime.voltage_exp, 'r-', 'MarkerSize', 8);
            hold on;
            plot(ftime.time_vec, V_ecm, 'b-', 'LineWidth', 1.5);        
            grid on;
            legend('Experimental / P2D', 'Optimized ECM', 'Location', 'best');
            title('Cell Voltage Fitting');
            xlabel('Time (s)');
            ylabel('Voltage (V)'); 
            
            % Subplot 3: Instantaneous Residual Error (V_exp - V_ecm)
            subplot(3,1,3);
            error_vec = ftime.voltage_exp - V_ecm;
            plot(ftime.time_vec, error_vec, 'g', 'LineWidth', 1.2);
            grid on;
            title('Residual Error (V_{exp} - V_{ecm})');
            xlabel('Time (s)');
            ylabel('Error (V)');
            
            % Text box displaying the global fitting score
            text_error = sprintf('Fitting Error (SSR): %.2e V²', fitting_error);
            text(0.05, 0.15, text_error, 'Units', 'normalized', ...
                 'BackgroundColor', 'white', ...   
                 'EdgeColor', 'black', ...        
                 'FontSize', 11, ...               
                 'FontWeight', 'bold');
        end

        function printResults(ftime, best_params, fitting_error)
            % Extraction des paramètres
            R0 = best_params(1);
            R1 = best_params(2);
            C1 = best_params(3);
            R2 = best_params(4);
            C2 = best_params(5);
            
            % Calcul des constantes de temps (tau = R * C)
            tau1 = R1 * C1;
            tau2 = R2 * C2;
            
            % Calcul de la RMSE pour un affichage plus parlant qu'une simple somme quadratique
            rmse_mv = sqrt(fitting_error / length(ftime.time_vec)) * 1000;
        
            fprintf('\n=======================================\n');
            fprintf('         FITTING RESULTS (TIME DOMAIN)     \n');
            fprintf('=======================================\n');
            fprintf('Total Error (SSR) : %.4e V²\n', fitting_error);

            fprintf('R0 : %.4e Ohms\n', R0);
            fprintf('R1 : %.4e Ohms\n', R1);
            fprintf('C1 : %.4e Farads\n', C1);
            fprintf('R2 : %.4e Ohms\n', R2);
            fprintf('C2 : %.4e Farads\n', C2);
        end
    end

end