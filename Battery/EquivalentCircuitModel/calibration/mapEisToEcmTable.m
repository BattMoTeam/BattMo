function ecm_table = mapEisToEcmTable()
    %% 1. Initialisation
    jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
    jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));
    
    jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                                   jsonstruct_geometry});
   
    
    [model, inputparams] = setupModelFromJson(jsonstruct);
    state_current = setupInitialState(model);
    
    soc_list = (1.0:-0.1:0.1)';
    N_steps = length(soc_list);
    
    results_matrix = zeros(N_steps, 6);
    
    freq = logspace(-4, 2, 50)./(2*pi);
    C_rate = 1;
    t_5percent = 3600 * 0.05 / C_rate;
    
    %% 2. Main loop
    for i = 1:N_steps
        current_soc = soc_list(i);
        fprintf('\n==================================================\n');
        fprintf('  EIS fitting for SOC = %.2f (%d/%d)\n', current_soc, i, N_steps);
        fprintf('==================================================\n');
        
        %% extraction de l'EIS par computeImpedance()
        options = struct();

        options.stateInitialization.computeSteadyState = false;

        options.stateInitialization.initializationSetup = 'given state'; 
        extrastructs.initstate = state_current;
        
        impsolv = ImpedanceSolver(inputparams, options, extrastructs);
        
        fprintf('Linearisation...\n');
        Z = impsolv.computeImpedance(freq);
        

        Z_re = real(Z)';
        Z_im = imag(Z)';
        
        omega = (2*pi).*freq;

        %% extraction des paramètres ECM
        p = runClassFitting(Z_re, Z_im, omega); 
        
        % Stockage dans la matrice de résultats
        results_matrix(i, 1) = current_soc;
        results_matrix(i, 2:6) = p;
        
        %% Décharge de 5% pour passer au SOC suivant
        if i < N_steps
            fprintf(' 5%% discharge...\n');
            inputparams.Control.controlPolicy = 'CCDischarge';
            inputparams.Control.DRate = C_rate;
            jsonstruct.TimeStepping.tmax = t_5percent;
            
            model.Control = model.setupControl(inputparams.Control);
            
           % CONVERSIONS MANUELLES DU SCHEDULE (Plus de setupScheduleStep !)
            % On découpe les 180s en 10 étapes de 18s pour aider le solveur à converger
            N_substeps = 10; 
            schedule_pulse = struct();
            schedule_pulse.step.val = ones(N_substeps, 1) * (t_5percent / N_substeps);
            schedule_pulse.step.control = ones(N_substeps, 1);
            % Cette fonction anonyme renvoie le courant Imax calculé par BattMo
            schedule_pulse.control.src = @(t, varargin) varargin{end}; 
            
            schedule_pulse.control.Control = struct('stopFunction', @(model, state, state0) false);

            % Simulation du pulse de décharge
            % simulateScheduleAD avec l'état en premier argument
            [~, states_pulse, ~] = simulateScheduleAD(state_current, model, schedule_pulse);

            % --- Micro-relaxation (Pause de 1 minute) ---
            fprintf(' 1 min pause...\n');
            inputparams.Control.DRate = 0; % Courant nul
            model.Control = model.setupControl(inputparams.Control);
            
            % Construction manuelle du schedule de repos
            N_rest_steps = 5;
            schedule_rest = struct();
            schedule_rest.step.val = ones(N_rest_steps, 1) * (60 / N_rest_steps);
            schedule_rest.step.control = ones(N_rest_steps, 1);
            schedule_rest.control.src = @(t, varargin) varargin{end};

            schedule_rest.control.Control = struct('stopFunction', @(model, state, state0) false);
            
            % Simulation du repos
            [~, states_rest, ~] = simulateScheduleAD(states_pulse{end}, model, schedule_rest);

            % L'état courant est mis à jour pour le prochain tour de boucle
            state_current = states_rest{end};
        end
    end
    
    %% 3. Printing final table
    ecm_table = array2table(results_matrix, 'VariableNames', ...
        {'SOC', 'R0', 'R1', 'C1', 'R2', 'C2'});
    
    fprintf('\n--- CARTOGRAPHIE DES PARAMÈTRES ECM TERMINÉE ---\n');
    disp(ecm_table);
    
    % automatic save
    save('ecm_map_results.mat', 'ecm_table');
end