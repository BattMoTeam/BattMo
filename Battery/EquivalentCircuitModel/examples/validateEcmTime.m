function validateEcmTime(file_name)
    
    %% ====================================================================
    %% STEP 1 : DATA LOADING AND OCV CATCHING
    %% ====================================================================
    
    % Set default file name if no argument is provided
    if nargin < 1
        file_name = 'ecm_map_results_10.mat'; 
    end

    % Check if the file exists and load it
    if exist(file_name, 'file')
        load(file_name, 'ecm_table');
    else
        error(['File not found: ', file_name]);
    end
    
    % Cathode (Positive)
    json_pe_path = fullfile('ParameterData','ParameterSets','Chen2020','chen2020_positive_electrode_interface.json');
    json_pe = parseBattmoJson(json_pe_path);
    [fn_ocp_pe, ~] = setupFunction(json_pe.openCircuitPotential);
    
    % Anode (Negative)
    json_ne_path = fullfile('ParameterData','ParameterSets','Chen2020','chen2020_negative_electrode_interface.json');
    json_ne = parseBattmoJson(json_ne_path);
    [fn_ocp_ne, ~] = setupFunction(json_ne.openCircuitPotential);
    
    % 1.3. Generation of the "Ground Truth": Continuous P2D discharge (BattMo)
    jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
    jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));
    jsonstruct = mergeJsonStructs({jsonstruct_material, jsonstruct_geometry});
    
    [model, inputparams] = setupModelFromJson(jsonstruct);
    
    %% New time-varying current profile setup
    custom_time    = [0, 100,  600, 800, 1500, 2000]; 
    custom_current = [0,   5,   -2,   0,    6,    0]; % Current profile in Amperes
    
    t_sim = custom_time(end); 
    
    % FIX HERE: Keeping CCDischarge ensures the MRST time-stepper steps correctly
    inputparams.Control.controlPolicy = 'CCDischarge'; 
    model.Control = model.setupControl(inputparams.Control);
    
    % Safe initialization using the standard policy branch
    state_init = setupInitialState(model);
    
    % Building the fine-grained time schedule (2-second steps)
    dt_step = 50; 
    N_sim_steps = floor(t_sim / dt_step);
    schedule_p2d = struct();
    schedule_p2d.step.val = ones(N_sim_steps, 1) * dt_step;
    schedule_p2d.step.control = ones(N_sim_steps, 1);
    
    % Dynamically feed the current profile into the simulation via interpolation
    schedule_p2d.control.src = @(t, varargin) interp1(custom_time, custom_current, t, 'previous', 0);
    schedule_p2d.control.Control = struct('stopFunction', @(model, state, state0) false);
    
    fprintf('Simulation of the P2D physical model (BattMo) in progress...\n');
    [~, states_p2d, ~] = simulateScheduleAD(state_init, model, schedule_p2d);
    
    % Extraction of time vectors from the P2D model
    N_points = length(states_p2d);
    time_vec = zeros(N_points, 1);
    V_p2d = zeros(N_points, 1);
    I_p2d = zeros(N_points, 1);
    soc_p2d = zeros(N_points, 1); % Pre-allocation du vrai SOC


    pe_guestStoichiometry0   = model.PositiveElectrode.Coating.ActiveMaterial.Interface.guestStoichiometry0;
    pe_guestStoichiometry100 = model.PositiveElectrode.Coating.ActiveMaterial.Interface.guestStoichiometry100;
    ne_guestStoichiometry0   = model.NegativeElectrode.Coating.ActiveMaterial.Interface.guestStoichiometry0;
    ne_guestStoichiometry100 = model.NegativeElectrode.Coating.ActiveMaterial.Interface.guestStoichiometry100;
    sat_conc_ne              = model.NegativeElectrode.Coating.ActiveMaterial.Interface.saturationConcentration;   
    
    for k = 1:N_points
        time_vec(k) = states_p2d{k}.time;
        V_p2d(k) = states_p2d{k}.Control.E;
        I_p2d(k) = states_p2d{k}.Control.I;
    end
    
    % Coulomb Counting to obtain SOC
     [Q_nominal, ~] = computeCellCapacity(model); 

    % Genuine SOC calculation tracking the real capacity depletion
    soc_p2d = 1 - cumtrapz(time_vec, I_p2d) / Q_nominal;

    % Security boundary to avoid interpolation out-of-bounds [0, 1]
    soc_p2d = max(0, min(1, soc_p2d));


    %% ====================================================================
    %% STEP 2 : ECM EQUATIONS IMPLEMENTATION (Simulation)
    %% ====================================================================
    fprintf('\n--- Step 2: Time-domain simulation of the ECM model (2-RC) ---\n');
    
    V_ecm = zeros(N_points, 1);
    
    % Initialization of ECM electrical state variables
    Uc1 = 0; % Initial voltage across RC branch #1 (V)
    Uc2 = 0; % Initial voltage across RC branch #2 (V)
    
    for k = 1:N_points
        if k == 1
            dt = time_vec(1);
        else
            dt = time_vec(k) - time_vec(k-1);
        end
        
        I = I_p2d(k);           
        soc = soc_p2d(k);       
        
        % 2.1. Dynamic interpolation of fitted ECM parameters according to SOC
        R0 = interp1(ecm_table.SOC, ecm_table.R0, soc, 'linear', 'extrap');
        R1 = interp1(ecm_table.SOC, ecm_table.R1, soc, 'linear', 'extrap');
        C1 = interp1(ecm_table.SOC, ecm_table.C1, soc, 'linear', 'extrap');
        R2 = interp1(ecm_table.SOC, ecm_table.R2, soc, 'linear', 'extrap');
        C2 = interp1(ecm_table.SOC, ecm_table.C2, soc, 'linear', 'extrap');
        
        
        x_pe = soc * (pe_guestStoichiometry100 - pe_guestStoichiometry0) + pe_guestStoichiometry0;              
        x_ne = soc * (ne_guestStoichiometry100 - ne_guestStoichiometry0) + ne_guestStoichiometry0;
        
        V_ocv = fn_ocp_pe(x_pe) - fn_ocp_ne(x_ne);
       
        % 2.3. Solving the Differential Equations of the RC branches (Explicit Euler)
        Uc1_next = Uc1 + (-Uc1 / (R1 * C1) + I / C1) * dt;
        Uc2_next = Uc2 + (-Uc2 / (R2 * C2) + I / C2) * dt;
        
        % 2.4. Calculation of the final cell voltage (Kirchhoff's Voltage Law)
        V_ecm(k) = V_ocv - R0 * I - Uc1 - Uc2;
        
        Uc1 = Uc1_next;
        Uc2 = Uc2_next;
    end
    
   %% ====================================================================
    %% STEP 3 : PLOTTING THE THREE FUNCTIONS (Visualization)
    %% ====================================================================
    fprintf('\n--- Step 3: Generating comparison plots ---\n');
    
    % Hauteur augmentée à 800 pixels pour accueillir confortablement les 3 subplots
    figure('Color', 'w', 'Name', 'Validation P2D vs ECM', 'Position', [100, 100, 900, 800]);
    
    % --- SUBPLOT 1: Tension (Axe gauche) & Courant (Axe droit) ---
    subplot(3,1,1);
    yyaxis left
    plot(time_vec, V_p2d, 'b-', 'LineWidth', 2.5); hold on;
    plot(time_vec, V_ecm, 'r--', 'LineWidth', 2);
    ylabel('Cell Voltage (V)', 'FontSize', 11);
    ax = gca; ax.YColor = 'k'; 
    
    yyaxis right
    plot(time_vec, I_p2d, 'g-', 'LineWidth', 1.5);
    ylabel('Current (A)', 'FontSize', 11);
    ax.YColor = 'g'; 
    
    grid on;
    title('Time-domain validation: Physical P2D Model vs Reduced ECM', 'FontSize', 13);
    legend('Ground Truth: BattMo P2D', 'Reduced Model: ECM (2-RC)', 'Input Current Profile', 'Location', 'SouthWest');
    
    % --- SUBPLOT 2: Profil du SOC physique au cours du temps ---
    subplot(3,1,2);
    plot(time_vec, soc_p2d * 100, 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 2); % Tracé en pourcentage
    grid on;
    ylabel('State of Charge (%)', 'FontSize', 11);
    title('Physical State of Charge (SOC) profile extracted from BattMo', 'FontSize', 11);
    
    % --- SUBPLOT 3: Erreur absolue de tension ---
    subplot(3,1,3);
    error_voltage = abs(V_p2d - V_ecm) * 1000; 
    plot(time_vec, error_voltage, 'k-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)', 'FontSize', 11);
    ylabel('Absolute Error (mV)', 'FontSize', 11);
    title('Instantaneous error of the equivalent circuit model', 'FontSize', 11);
    
    rmse = sqrt(mean((V_p2d - V_ecm).^2)) * 1000;
    fprintf('*** Analysis complete! Root Mean Square Error (RMSE): %.2f mV ***\n', rmse);
end