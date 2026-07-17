function jsonstruct_ecm = calibrateEcmFromP2D(jsonstruct_p2d, soc_step)
    % CALIBRATEECMFROMP2D Calibrates a 2-RC Equivalent Circuit Model (ECM)
    % from a physics-based P2D model across a range of State of Charge (SOC) points.
    %
    % Inputs:
    %   jsonstruct_p2d - P2D model configuration structure
    %   soc_step       - SOC decrement interval (default: 0.1)
    %
    % Output:
    %   jsonstruct_ecm - Standardized BattMo ECM configuration structure
    
    %% 1. Initialization
    if nargin < 2
        soc_step = 0.1; 
    end
    
    [model, inputparams] = setupModelFromJson(jsonstruct_p2d);
    state_current = setupInitialState(model);
    
    soc_list = (1.0:-soc_step:0.1)';
    N_steps = length(soc_list);
    
    results_matrix = zeros(N_steps, 6); % Stores [SOC, R0, R1, C1, R2, C2]
    ocv_list = zeros(N_steps, 1);       % Stores Open Circuit Voltage for each SOC
    
    freq = logspace(-4, 2, 50)./(2*pi);
    C_rate = 1;
    t_xpercent = 3600 * soc_step / C_rate;
    
    %% 2. Main Loop
    for i = 1:N_steps
        current_soc = soc_list(i);
        fprintf('\n==================================================\n');
        fprintf('  EIS fitting for SOC = %.2f (%d/%d)\n', current_soc, i, N_steps);
        fprintf('==================================================\n');
        
        % Robustness fallback: ensure the initial state has evaluated terminal voltage (.Control.E)
        if i == 1 && (~isfield(state_current, 'Control') || ~isfield(state_current.Control, 'E'))
            inputparams.Control.DRate = 0;
            model.Control = model.setupControl(inputparams.Control);
            schedule_init = struct();
            schedule_init.step.val = 0.001;
            schedule_init.step.control = 1;
            schedule_init.control.src = @(t, varargin) 0;
            schedule_init.control.Control = struct('stopFunction', @(model, state, state0) false);
            [~, states_init, ~] = simulateScheduleAD(state_current, model, schedule_init);
            state_current = states_init{end};
        end
        
        % Capture the relaxed terminal voltage as the Open Circuit Voltage (OCV)
        ocv_list(i) = state_current.Control.E;
        fprintf('  Measured Open Circuit Voltage (OCV): %.4f V\n', ocv_list(i));
        
        % Run the Impedance Solver
        options = struct();
        options.stateInitialization.computeSteadyState = false;
        options.stateInitialization.initializationSetup = 'given state'; 
        extrastructs.initstate = state_current;
        
        impsolv = ImpedanceSolver(inputparams, options, extrastructs);
        
        fprintf('Linearization...\n');
        Z = impsolv.computeImpedance(freq);
        Z_re = real(Z)';
        Z_im = imag(Z)';
        
        omega = (2*pi).*freq;
        
        % Run parameter identification
        eisdata = struct('omega', omega, ...
                 'Z_re_exp', Z_re, ...
                 'Z_im_exp', Z_im);


        params0 = [0.05052; 1.12673; 59119.9; 0.03155; 11054.0];  % initial condition: C1>2*C2
        
        feis = FittingEIS(eisdata, params0);
        [p, ~] = feis.run();

        
        % Store values in the matrix
        results_matrix(i, 1) = current_soc;
        results_matrix(i, 2:6) = p;
        
        % Discharge step to reach the next SOC target
        if i < N_steps
            fprintf(' %02d %% discharge...\n', round(soc_step * 100));
            inputparams.Control.controlPolicy = 'CCDischarge';
            inputparams.Control.DRate = C_rate;
            jsonstruct_p2d.TimeStepping.tmax = t_xpercent;
            
            model.Control = model.setupControl(inputparams.Control);
            
            N_substeps = 10; 
            schedule_pulse = struct();
            schedule_pulse.step.val = ones(N_substeps, 1) * (t_xpercent / N_substeps);
            schedule_pulse.step.control = ones(N_substeps, 1);
            schedule_pulse.control.src = @(t, varargin) varargin{end}; 
            schedule_pulse.control.Control = struct('stopFunction', @(model, state, state0) false);
            
            [~, states_pulse, ~] = simulateScheduleAD(state_current, model, schedule_pulse);
            
            % Relaxation step (1 minute)
            fprintf(' 1 min pause...\n');
            inputparams.Control.DRate = 0; 
            model.Control = model.setupControl(inputparams.Control);
            
            N_rest_steps = 5;
            schedule_rest = struct();
            schedule_rest.step.val = ones(N_rest_steps, 1) * (60 / N_rest_steps);
            schedule_rest.step.control = ones(N_rest_steps, 1);
            schedule_rest.control.src = @(t, varargin) varargin{end};
            schedule_rest.control.Control = struct('stopFunction', @(model, state, state0) false);
            
            [~, states_rest, ~] = simulateScheduleAD(states_pulse{end}, model, schedule_rest);
            state_current = states_rest{end};
        end
    end
    
    %% 3. Post-Processing & Sorting (Ascending SOC Order)
    % BattMo ECM JSON requires OCP and parameters defined with ascending SOC (0 -> 1)
    [~, sort_idx] = sort(soc_list, 'ascend');
    soc_asc = soc_list(sort_idx)';
    ocv_asc = ocv_list(sort_idx)';
    
    results_matrix_sorted = results_matrix(sort_idx, :);
    R0_asc = results_matrix_sorted(:, 2)';
    R1_asc = results_matrix_sorted(:, 3)';
    C1_asc = results_matrix_sorted(:, 4)';
    R2_asc = results_matrix_sorted(:, 5)';
    C2_asc = results_matrix_sorted(:, 6)';
    
    %% 4. Build jsonstruct_ecm
    jsonstruct_ecm = struct();
    
    % Nominal capacity calculation (convert from Ampere-seconds to Ampere-hours)
    [Q_nominal, ~] = computeCellCapacity(model); 
    jsonstruct_ecm.nominalcellcapacity = Q_nominal / 3600; 
    
    jsonstruct_ecm.initSOC = 1.0;
    
    % Extract lower voltage cutoff if present in the P2D model, otherwise use default
    if isfield(jsonstruct_p2d.Control, 'lowerCutoffVoltage')
        jsonstruct_ecm.lowerVoltageCutoff = jsonstruct_p2d.Control.lowerCutoffVoltage;
    else
        jsonstruct_ecm.lowerVoltageCutoff = 3.0; 
    end
    
    jsonstruct_ecm.totalTime = 4200;
    
    % Tabulated input current profiles (aligned with template)
    jsonstruct_ecm.I.functionFormat = 'tabulated';
    jsonstruct_ecm.I.argumentList = {'time'};
    jsonstruct_ecm.I.dataX = [0, 100, 100.1, 120, 120.1, 300, 300.1, 3800, 3800.1, 4000, 4000.1, 4020, 4020.1, 4200];
    jsonstruct_ecm.I.dataY = [0, 0,   60,    60,  0,     0,   60,    60,   0,      0,    60,     60,   0,      0];
    
    % Tabulated Open Circuit Potential (OCP)
    jsonstruct_ecm.OCP.functionFormat = 'tabulated';
    jsonstruct_ecm.OCP.argumentList = {'SOC'};
    jsonstruct_ecm.OCP.dataX = soc_asc;
    jsonstruct_ecm.OCP.dataY = ocv_asc;
    
    % Tabulated R0
    jsonstruct_ecm.R0.functionFormat = 'tabulated';
    jsonstruct_ecm.R0.argumentList = {'SOC'};
    jsonstruct_ecm.R0.dataX = soc_asc;
    jsonstruct_ecm.R0.dataY = R0_asc;
    
    % Tabulated R1 & C1 (First RC branch)
    jsonstruct_ecm.R1.functionFormat = 'tabulated';
    jsonstruct_ecm.R1.argumentList = {'SOC'};
    jsonstruct_ecm.R1.dataX = soc_asc;
    jsonstruct_ecm.R1.dataY = R1_asc;
    
    jsonstruct_ecm.C1.functionFormat = 'tabulated';
    jsonstruct_ecm.C1.argumentList = {'SOC'};
    jsonstruct_ecm.C1.dataX = soc_asc;
    jsonstruct_ecm.C1.dataY = C1_asc;
    
    % Tabulated R2 & C2 (Second RC branch)
    jsonstruct_ecm.R2.functionFormat = 'tabulated';
    jsonstruct_ecm.R2.argumentList = {'SOC'};
    jsonstruct_ecm.R2.dataX = soc_asc;
    jsonstruct_ecm.R2.dataY = R2_asc;
    
    jsonstruct_ecm.C2.functionFormat = 'tabulated';
    jsonstruct_ecm.C2.argumentList = {'SOC'};
    jsonstruct_ecm.C2.dataX = soc_asc;
    jsonstruct_ecm.C2.dataY = C2_asc;
    
    % Print and save summary
    ecm_table = array2table(results_matrix_sorted, 'VariableNames', ...
        {'SOC', 'R0', 'R1', 'C1', 'R2', 'C2'});
    disp(ecm_table);
    
    file_name = sprintf('jsonstruct_ecm_calibrated_%02d.mat', round(soc_step * 100));
    save(file_name, 'jsonstruct_ecm');
    fprintf('jsonstruct_ecm structure saved successfully in: %s\n', file_name);

    % % Convert the struct to a formatted JSON string and display it
    % json_readable = jsonencode(jsonstruct_ecm, 'PrettyPrint', true);
    % disp(json_readable);
end