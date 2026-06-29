function ecm_table = mapEisToEcmTable()
    %% 1. Initialisation
    json_path = 'C:\Users\Alexandre Fichter\Documents\stage_3A\contenu stage\data_August\jp3-params\jp3-opt-1d-full.json';
    jsonstruct = parseBattmoJson(json_path);
    [model, inputparams] = setupModelFromJson(jsonstruct);
    state_current = setupInitialState(model);
    
    soc_list = (1.0:-0.05:0.05)';
    N_steps = length(soc_list);
    
    results_matrix = zeros(N_steps, 6);
    

    C_rate = 1;
    t_5percent = 3600 * 0.05 / C_rate;
    
    %% 2. Main loop
    for i = 1:N_steps
        current_soc = soc_list(i);
        fprintf('\n==================================================\n');
        fprintf('  EIS fitting for SOC = %.2f (%d/%d)\n', current_soc, i, N_steps);
        fprintf('==================================================\n');
        
        %% extraction de l'EIS par computeImpedance()
        
        %% extraction des paramètres ECM
        p = runClassFitting(Z_re, Z_im, omega); 
        
        % Stockage dans la matrice de résultats
        results_matrix(i, 1) = current_soc;
        results_matrix(i, 2:6) = p; % Insère [R0, R1, C1, R2, C2]
        
        %% Décharge de 5% pour passer au SOC suivant
        if i < N_steps
            fprintf(' 5%% discharge...\n');
            inputparams.Control.controlPolicy = 'CCDischarge';
            inputparams.Control.DRate = C_rate;
            jsonstruct.TimeStepping.tmax = t_5percent;
            
            schedule_pulse = model.setupSchedule(jsonstruct);
            [~, states_pulse] = simulateBattery(model, state_current, schedule_pulse);
            
            % Micro-relaxation 
            inputparams.Control.DRate = 0;
            jsonstruct.TimeStepping.tmax = 60; % 1 minute de repos
            schedule_rest = model.setupSchedule(jsonstruct);
            [~, states_rest] = simulateBattery(model, states_pulse{end}, schedule_rest);
            
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