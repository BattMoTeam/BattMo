classdef FittingTime

    % FittingTime - Fit a 2-RC Thevenin ECM to P2D simulated time-domain data
    %
    % Usage:
    %   ftime = FittingTime(jsonstruct, time_vec, current_exp, params0, scales)
    %
    % Inputs:
    %   jsonstruct  - BattMo JSON structure used to setup P2D model and OCV
    %   time_vec    - vector of time instants (s) for experimental current
    %   current_exp - experimental current vector (A) matched to time_vec
    %   voltage_exp - experimental voltage vector (V) matched to time_vec
    %   params0     - initial guess [R0; R1; C1; R2; C2]
    %   scales      - normalization bounds vector: [pmin(1:5); pmax(1:5)]
    %
    % Main methods:
    %   optimizationBFGS(ftime)       - run BFGS optimization in unit box
    %   optifunc(ftime, p)            - objective and gradient for params p
    %   plotresults_thevenin(...)     - plot fitted vs reference voltages
    %   printResults(...)             - print numeric results and diagnostics
    %   plottest(...)                 - validate ECM on a new current profile
    %
    % Notes:
    %   - The class simulates a BattMo P2D model to produce reference OCV,
    %     voltage and current traces which are used as the 'experimental'
    %     data for fitting the equivalent-circuit model parameters.
    %   - Parameter vector is [R0, R1, C1, R2, C2]. Scales specify lower and
    %     upper bounds used for normalization between 0 and 1.
    %
    % Example:
    %   ftime = FittingTime(jsonstruct, t, Iexp, params0, scales);
    %   [v, hist, best, err] = ftime.optimizationBFGS();
    %   ftime.plotresults_thevenin(best, err);



    properties(SetAccess = immutable)

        inputs
        
    end
    
    properties

        %% experimental value
        time_vec     % time steps for current_exp and voltage_exp
        current_exp  % experimental current values
        voltage_exp  % experimental voltage values
        initialSOC   % soc at start of experiment
        capacity     % battery capacity, needed to compute soc
        
        %% OCV values
        ocvFunction % open circuit voltage as a function of soc

        %% P2D model parameters
        jsonstruct  % input structure for P2D model which gives us the the ocv and the experimental data
        
        %% optimization parameters
        params0     % initial guess for the parameters to be fitted
        scales      % scaling for the unit box optimization

        %%
        options

    end

    methods        
    
        function ftime = FittingTime(params0, inputs, options)

            if nargin < 3
                options = [];
            end

            ftime.inputs = inputs;

            options = setDefaultStructField(options, 'useDefaultScales', true);
            options = setDefaultStructField(options, 'useP2Dmodel', false);
            options = setDefaultStructField(options, 'scaleFactor', 1000);
            options = setDefaultStructField(options, 'useSimulationOutput', false);
            
            ftime.options = options;

            ftime.params0 = params0;
            
            if ftime.options.useP2Dmodel
                ftime.jsonstruct = inputs.jsonstruct;
                ftime = ftime.setupFromP2D();
            else
                ftime.time_vec    = inputs.time_vec;
                ftime.current_exp = inputs.current_exp;
                ftime.voltage_exp = inputs.voltage_exp;
                ftime.initialSOC  = inputs.initialSOC;
                ftime.ocvFunction = inputs.ocvFunction;
            end

            if ftime.options.useDefaultScales
                ftime = ftime.setupDefaultScales(options.scaleFactor);
            else
                ftime.scales = inputs.scales;
            end
            

        end

        function ftime = setupDefaultScales(ftime, factor)

            params0 = ftime.params0;
            
            pmin = params0 / factor;
            pmax = params0 * factor;
            scales = [pmin; pmax];

            ftime.scales = scales;
            
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

        function ftime = setupFromP2D(ftime)

            jsonstruct  = ftime.jsonstruct;

            if ftime.options.useSimulationOutput
                output = ftime.inputs.output;
            else
                output = runBattery(jsonstruct);
            end
            
            model = output.model;

            ftime.time_vec      = cellfun(@(state) state.time, output.states);
            ftime.current_exp   = cellfun(@(state) state.Control.I, output.states);
            ftime.voltage_exp   = cellfun(@(state) state.Control.E, output.states);

            ftime.initialSOC = model.SOC;
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            am  = 'ActiveMaterial';
            itf = 'Interface';
            
            ne_ocp = model.(ne).(co).(am).(itf).computeOCP;
            ne_100 = model.(ne).(co).(am).(itf).guestStoichiometry100;
            ne_0   = model.(ne).(co).(am).(itf).guestStoichiometry0;
            pe_ocp = model.(pe).(co).(am).(itf).computeOCP;
            pe_100 = model.(pe).(co).(am).(itf).guestStoichiometry100;
            pe_0   = model.(pe).(co).(am).(itf).guestStoichiometry0;
            
            ftime.ocvFunction = @(soc) (pe_ocp(soc*(pe_100 - pe_0) + pe_0) - ...
                                        ne_ocp(soc* (ne_100 - ne_0) + ne_0));
            
            ftime.capacity = computeCellCapacity(model);
            
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

            [min_value, best_params_norm, history] = unitBoxBFGS(params0_norm               , ...                             
                                                                 f_opt                      , ...           
                                                                 'maximize'       , false   , ...      
                                                                 'linIneq'        , struct('A', A_custom, 'b', b_custom),  ...   %A*u<=b
                                                                 'enforceFeasible', true    , ...     
                                                                 'maxIt'          , max_iter, ...               
                                                                 'objChangeTol'   , tol_obj , ...       
                                                                 'gradTol'        , tol_grad, ...
                                                                 'lineSearchMaxIt', 100);
            
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

            soc = ftime.initialSOC;
            Q   = ftime.capacity;
            
            N_points = length(ftime.time_vec);

            V_ecm = zeros(N_points, 1);

            Uc1 = 0; 
            Uc2 = 0;
            
            % Sensibilities initialisation 

            % we do not need s_R0 as it is trivial
            s_R1 = 0; % dUc1 / dR1
            s_C1 = 0; % dUc1 / dC1
            s_R2 = 0; % dUc2 / dR2
            s_C2 = 0; % dUc2 / dC2
            
            g_true = zeros(5, 1);
            
            for k = 1 : N_points
                
                if k == 1
                    dt = ftime.time_vec(1);
                else
                    dt = ftime.time_vec(k) - ftime.time_vec(k-1);
                end
                
                I     = ftime.current_exp(k);           
                V_ocv = ftime.ocvFunction(soc); 
                
                % 1. Tension du modèle au pas k
                V_ecm(k) = V_ocv - R0 * I - Uc1 - Uc2;
                
                % 2. Erreur instantanée
                err = ftime.voltage_exp(k) - V_ecm(k);

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
                soc = soc - I*dt/Q;

            end
            
            % Valeur de la fonction objectif (Erreur quadratique totale)
            v = sum((ftime.voltage_exp - V_ecm).^2);

            if isnan(v) || isinf(v) || (nargout > 1 && (any(isnan(g_true)) || any(isinf(g_true))))
                v = 1e10; % Assigne un coût massif pour rejeter le point instable
                if nargout > 1
                    g_norm = zeros(5, 1); % Renvoie un gradient plat pour forcer le pas arrière
                end
                return;
            end

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

            soc = ftime.initialSOC;
            Q   = ftime.capacity;
            
            for k = 1 : N_points
                
                if k == 1
                    dt = ftime.time_vec(1);
                else
                    dt = ftime.time_vec(k) - ftime.time_vec(k-1);
                end
                I = ftime.current_exp(k);
                
                V_ocv = ftime.ocvFunction(soc);
                
                V_ecm(k) = V_ocv - R0 * I - Uc1 - Uc2;
                
                Uc1 = Uc1 + (-Uc1 / (R1 * C1) + I / C1) * dt;
                Uc2 = Uc2 + (-Uc2 / (R2 * C2) + I / C2) * dt;

                soc = soc - I*dt/Q;
                
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


        function plottest(ftime, best_params, current_test)

            % utility function to compare with a given modified current input, see definition below
            
            % 1. Recalculate the final voltage trajectory using best_params
            p_safe = max(best_params, 1e-8);
            
            R0 = p_safe(1);
            R1 = p_safe(2);
            C1 = p_safe(3);
            R2 = p_safe(4);
            C2 = p_safe(5);
            
            % Get original training time grid as a column vector
            time_test = ftime.time_vec(:); 
            current_test = current_test(:);
            
            [voltage_test, ocv_test, time_test_sim, current_test_sim] = ftime.setupOCP(time_test, current_test);
    
            % Force all outputs to be column vectors to avoid any plot orientation errors
            voltage_test     = voltage_test(:);
            ocv_test         = ocv_test(:);
            time_test_sim    = time_test_sim(:);
            current_test_sim = current_test_sim(:);
            
            % Synchronize the simulation length with the actual P2D points (e.g., 99)
            N_points = length(time_test_sim);
            
            % Initialize ECM state variables
            V_ecm = zeros(N_points, 1);
            Uc1 = 0; 
            Uc2 = 0;

            for k = 1:N_points
                if k == 1
                    dt = time_test_sim(1);
                else
                    dt = time_test_sim(k) - time_test_sim(k-1);
                end
                I = current_test_sim(k); % Courant réel synchronisé simulé par BattMo
                V_ocv = ocv_test(k);     % OCV réel synchronisé simulé par BattMo
                V_ecm(k) = V_ocv - R0 * I - Uc1 - Uc2;
                
                Uc1 = Uc1 + (-Uc1 / (R1 * C1) + I / C1) * dt;
                Uc2 = Uc2 + (-Uc2 / (R2 * C2) + I / C2) * dt;
            end
            
            % 2. Generate time-domain plots using the synchronized time_test grid
            figure('Name', 'ECM Time Domain Validation Results', 'NumberTitle', 'off');
            
            % Subplot 1: Applied Current Profile
            subplot(3,1,1);
            plot(time_test_sim, current_test_sim, 'k', 'LineWidth', 1.5);
            grid on;
            title('Applied Current Profile (P2D Simulated Range)');
            xlabel('Time (s)');
            ylabel('Current (A)');
            
            % Subplot 2: Voltage Comparison (P2D vs ECM)
            subplot(3,1,2);
            plot(time_test_sim, voltage_test, 'r-', 'LineWidth', 1.5); 
            hold on;
            plot(time_test_sim, V_ecm, 'b--', 'LineWidth', 1.5);        
            grid on;
            legend('Ground Truth: P2D Model', 'Validated ECM', 'Location', 'best');
            title('Cell Voltage Validation');
            xlabel('Time (s)');
            ylabel('Voltage (V)'); 
            
            % Subplot 3: Instantaneous Residual Error (V_p2d - V_ecm)
            subplot(3,1,3);
            error_vec = voltage_test - V_ecm; 
            plot(time_test_sim, error_vec, 'g', 'LineWidth', 1.2);
            grid on;
            title('Residual Error (V_{p2d} - V_{ecm})');
            xlabel('Time (s)');
            ylabel('Error (V)');
            
        end
    end

end
