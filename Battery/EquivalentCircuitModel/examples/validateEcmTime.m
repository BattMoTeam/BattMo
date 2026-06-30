function validateEcmTime()
    %% ====================================================================
    %% STEP 1 : DATA LOADING AND OCV CATCHING
    %% ====================================================================
    
    if exist('ecm_map_results.mat', 'file')
        load('ecm_map_results.mat', 'ecm_table');
    else
        error('file not found ecm_map_results.mat. Run mapEisToEcmTable first.');
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
    % We run a constant current discharge WITHOUT any pauses to serve as a temporal reference
    jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
    jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));
    jsonstruct = mergeJsonStructs({jsonstruct_material, jsonstruct_geometry});
    
    [model, inputparams] = setupModelFromJson(jsonstruct);
    state_init = setupInitialState(model);
    
    C_rate = 1; 
    t_sim = 3200; % Duration in seconds, be careful not to reach the computed SOC limit (here <10%)
    
    inputparams.Control.controlPolicy = 'CCDischarge';
    inputparams.Control.DRate = C_rate;
    model.Control = model.setupControl(inputparams.Control);
    
    % Continuous schedule construction
    dt_step = 10;
    N_sim_steps = floor(t_sim / dt_step);
    schedule_p2d = struct();
    schedule_p2d.step.val = ones(N_sim_steps, 1) * dt_step;
    schedule_p2d.step.control = ones(N_sim_steps, 1);
    schedule_p2d.control.src = @(t, varargin) varargin{end};
    schedule_p2d.control.Control = struct('stopFunction', @(model, state, state0) false);
    
    fprintf('Simulation of the P2D physical model (BattMo) in progress...\n');
    [~, states_p2d, ~] = simulateScheduleAD(state_init, model, schedule_p2d);
    
    % Extraction of time vectors from the P2D model
    N_points = length(states_p2d);
    time_vec = zeros(N_points, 1);
    V_p2d = zeros(N_points, 1);
    I_p2d = zeros(N_points, 1);
    
    for k = 1:N_points
        time_vec(k) = states_p2d{k}.time;
        V_p2d(k) = states_p2d{k}.Control.E;
        I_p2d(k) = states_p2d{k}.Control.I;
    end
    
    % Coulomb Counting to obtain SOC
    total_charge = trapz(time_vec, I_p2d);
    soc_p2d = 1 - cumtrapz(time_vec, I_p2d) / total_charge;
    
    %% ====================================================================
    %% STEP 2 : ECM EQUATIONS IMPLEMENTATION (Simulation)
    %% ====================================================================
    fprintf('\n--- Step 2: Time-domain simulation of the ECM model (2-RC) ---\n');
    
    V_ecm = zeros(N_points, 1);
    
    % Initialization of ECM electrical state variables
    Uc1 = 0; % Initial voltage across RC branch #1 (V)
    Uc2 = 0; % Initial voltage across RC branch #2 (V)
    
    for k = 1:N_points
        % Calculation of the real time-step (dt)
        if k == 1
            dt = time_vec(1);
        else
            dt = time_vec(k) - time_vec(k-1);
        end
        
        I = I_p2d(k);           % Current measured at this timestamp (A)
        soc = soc_p2d(k);       % SOC estimated at this timestamp
        
        % 2.1. Dynamic interpolation of fitted ECM parameters according to SOC
        R0 = interp1(ecm_table.SOC, ecm_table.R0, soc, 'linear', 'extrap');
        R1 = interp1(ecm_table.SOC, ecm_table.R1, soc, 'linear', 'extrap');
        C1 = interp1(ecm_table.SOC, ecm_table.C1, soc, 'linear', 'extrap');
        R2 = interp1(ecm_table.SOC, ecm_table.R2, soc, 'linear', 'extrap');
        C2 = interp1(ecm_table.SOC, ecm_table.C2, soc, 'linear', 'extrap');
        
        % 2.2. Evaluation of the exact thermodynamic OCV via setupFunction
        % For the Positive Electrode (NMC Cathode)
        pe_guestStoichiometry0   = model.PositiveElectrode.Coating.ActiveMaterial.Interface.guestStoichiometry0;
        pe_guestStoichiometry100 = model.PositiveElectrode.Coating.ActiveMaterial.Interface.guestStoichiometry100;
        
        % For the Negative Electrode (Graphite Anode)
        ne_guestStoichiometry0   = model.NegativeElectrode.Coating.ActiveMaterial.Interface.guestStoichiometry0;
        ne_guestStoichiometry100 = model.NegativeElectrode.Coating.ActiveMaterial.Interface.guestStoichiometry100;
        
        % Note: Chen2020 uses stoichiometry (0 to 1) instead of raw SOC.
        x_pe = soc * (pe_guestStoichiometry100 - pe_guestStoichiometry0) + pe_guestStoichiometry0;              
        % Translation of SOC into stoichiometry for the anode (Negative)
        x_ne = soc * (ne_guestStoichiometry100 - ne_guestStoichiometry0) + ne_guestStoichiometry0;
        
        % 3. Calculation of the true open circuit voltage of the cell
        V_ocv = fn_ocp_pe(x_pe) - fn_ocp_ne(x_ne);
       
        % 2.3. Solving the Differential Equations of the RC branches (Explicit Euler)
        Uc1_next = Uc1 + (-Uc1 / (R1 * C1) + I / C1) * dt;
        Uc2_next = Uc2 + (-Uc2 / (R2 * C2) + I / C2) * dt;
        
        % 2.4. Calculation of the final cell voltage (Kirchhoff's Voltage Law)
        V_ecm(k) = V_ocv - R0 * I - Uc1 - Uc2;
        
        % Passing states to the next step
        Uc1 = Uc1_next;
        Uc2 = Uc2_next;
    end
    
    %% ====================================================================
    %% STEP 3 : PLOTTING THE TWO FUNCTIONS (Visualization)
    %% ====================================================================
    fprintf('\n--- Step 3: Generating comparison plots ---\n');
    
    figure('Color', 'w', 'Name', 'Validation P2D vs ECM', 'Position', [100, 100, 900, 600]);
    
    % Voltage Plot
    subplot(2,1,1);
    plot(time_vec, V_p2d, 'b-', 'LineWidth', 2.5); hold on;
    plot(time_vec, V_ecm, 'r--', 'LineWidth', 2);
    grid on;
    ylabel('Cell Voltage (V)', 'FontSize', 11);
    title('Time-domain validation: Physical P2D Model vs Reduced ECM', 'FontSize', 13);
    legend('Ground Truth : BattMo P2D', 'Reduced Model : ECM (2-RC fitted by EIS)', 'Location', 'SouthWest');
    
    % Absolute Voltage Error Plot
    subplot(2,1,2);
    error_voltage = abs(V_p2d - V_ecm) * 1000; % Conversion to mV
    plot(time_vec, error_voltage, 'k-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)', 'FontSize', 11);
    ylabel('Absolute Error (mV)', 'FontSize', 11);
    title('Instantaneous error of the equivalent circuit model', 'FontSize', 11);
    
    % Calculation of the global mean error for the analysis report
    rmse = sqrt(mean((V_p2d - V_ecm).^2)) * 1000;
    fprintf('*** Analysis complete! Root Mean Square Error (RMSE): %.2f mV ***\n', rmse);
end