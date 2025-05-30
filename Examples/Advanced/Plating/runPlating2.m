%% run stand-alone active material model with lithium plating for different kInter values

clear
close all

%% Define kInter values to test
kInter_values = [6e3, 6e4, 6e5];
n_simulations = length(kInter_values);

% Initialize storage for results
results = cell(n_simulations, 1);
all_states = cell(n_simulations, 1);
all_time = cell(n_simulations, 1);

%% Loop over different kInter values
for sim_idx = 1:n_simulations
    
    fprintf('Running simulation %d/%d with kInter = %.0e\n', sim_idx, n_simulations, kInter_values(sim_idx));
    
    %% Setup the properties of Li-ion battery materials and cell design
    
    jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
    
    inputparams = BatteryInputParams(jsonstruct);
    
    use_cccv = true;
    if use_cccv
        inputparams.SOC = 0;
        cccvstruct = struct( 'controlPolicy'     , 'CCCV'       , ...
                             'initialControl'    , 'charging', ...
                             'numberOfCycles'    , 2            , ...
                             'CRate'             , 1.5          , ...
                             'DRate'             , 1            , ...
                             'lowerCutoffVoltage', 3            , ...
                             'upperCutoffVoltage', 4            , ...
                             'dIdtLimit'         , 1e-2         , ...
                             'dEdtLimit'         , 1e-4);
        cccvinputparams = CcCvControlModelInputParams(cccvstruct);
        inputparams.Control = cccvinputparams;
    end
    
    ne      = 'NegativeElectrode';
    pe      = 'PositiveElectrode';
    elyte   = 'Electrolyte';
    thermal = 'ThermalModel';
    co      = 'Coating';
    am      = 'ActiveMaterial';
    itf     = 'Interface';
    sd      = 'SolidDiffusion';
    ctrl    = 'Control';
    cc      = 'CurrentCollector';
    
    jsonstruct.use_thermal = false;
    jsonstruct.include_current_collectors = false;
    
    jsonstruct.(ne).(co).(am).diffusionModelType = 'full';
    jsonstruct.(pe).(co).(am).diffusionModelType = 'full';
    
    jsonstruct.(ne).(co).(am).useLithiumPlating = true;
    
    % Flag pour modèle stand-alone
    jsonstruct.(ne).(co).(am).isRootSimulationModel = true;
    
    jsonstruct_lithium_plating = parseBattmoJson(fullfile('Examples', 'Advanced', 'Plating', 'lithium_plating.json'));
    
    % Modify kInter value for current simulation
    jsonstruct_lithium_plating.LithiumPlating.kInter = kInter_values(sim_idx);
    
    jsonstruct.(ne).(co).(am).LithiumPlating = jsonstruct_lithium_plating.LithiumPlating;
    
    scenario = 'charge';
    
    %% following is not used at particle level (but necessary to initiliaze full battery below)
    switch scenario
      case 'charge'
        jsonstruct.Control.controlPolicy = 'CCCharge';
      case 'discharge'
        jsonstruct.Control.controlPolicy = 'CCDischarge';
      otherwise
        error('scenario not recognized');
    end
    
    % Setup InputParams
    inputparams = BatteryInputParams(jsonstruct);
    inputparams = inputparams.(ne).(co).(am);
    inputparams.Interface.openCircuitPotential.functionname = 'computeOCP_Graphite_Latz';
    inputparams.LithiumPlating.volumeFraction = inputparams.SolidDiffusion.volumeFraction;
    
    %% Setup the model
    
    model = ActiveMaterial(inputparams);
    
    %% Equip model for simulation
    model = model.setupForSimulation();
    
    %% Setup initial state
    
    sd  = 'SolidDiffusion';
    itf = 'Interface';
    
    cElectrolyte   = 5e-1*mol/litre;
    phiElectrolyte = 0;
    T              = 298;
    
    switch scenario
      case 'charge'
        cElectrodeInit = (model.(itf).guestStoichiometry0)*(model.(itf).saturationConcentration);
      case 'discharge'
        cElectrodeInit = (model.(itf).guestStoichiometry100)*(model.(itf).saturationConcentration);
      otherwise
        error('scenario not recognized');
    end
    
    N = model.(sd).N;
    initState.(sd).c        = cElectrodeInit*ones(N, 1);
    initState.(sd).cSurface = cElectrodeInit;
    
    initState.T = T;
    initState.(itf).cElectrolyte   = cElectrolyte;
    initState.(itf).phiElectrolyte = phiElectrolyte;
    
    initState = model.evalVarName(initState, {itf, 'OCP'});
    OCP = initState.(itf).OCP;
    initState.E = OCP + phiElectrolyte;
    
    lp = 'LithiumPlating';
    F = model.LithiumPlating.F;
    R = model.LithiumPlating.R;
    
    nPl0 = model.LithiumPlating.nPl0;
    r = model.LithiumPlating.particleRadius;
    poros = model.LithiumPlating.volumeFraction;
    platedConcentration0 = nPl0 * poros / ((4/3)*pi*r^3);
    
    initState.(lp).platedConcentration = platedConcentration0/(exp((F*OCP)/(R*T)) - 1)^(1/4);
    initState.(lp).phiSolid       = initState.E;
    initState.(lp).phiElectrolyte = phiElectrolyte;
    initState.(lp).cElectrolyte   = cElectrolyte;
    initState.(lp).nSEI = 0;
    
    %% setup schedule
    
    Iref = 5e-12; % calibrated set to work on this example
    Imax = 5e1*Iref;
    total = 1*hour*(Iref/Imax);
    n     = 100;
    dt    = total/n;
    step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));
    
    tup = 1*second*(Iref/Imax);
    
    switch scenario
      case 'charge'
        srcfunc = @(time) rampupControl(time, tup, -Imax); %0 pour tourner à vide
      case 'discharge'
        srcfunc = @(time) rampupControl(time, tup, Imax); %0 pour tourner à vide
      otherwise
        error('scenario not recognized');
    end
    
    cmin = (model.(itf).guestStoichiometry0)*(model.(itf).saturationConcentration);
    control.stopFunction = @(model, state, state0_inner) (state.(sd).cSurface <= cmin);
    control.src = srcfunc;
    
    schedule = struct('control', control, 'step', step);
    
    %% setup non-linear solver
    
    nls = NonLinearSolver();
    nls.errorOnFailure = false;
    model.nonlinearTolerance = 1e-2;
    
    %% Run simulation
    
    model.verbose = false; % Reduce verbosity for multiple simulations
    [~, states, report] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);
    
    %% Process results
    
    ind = cellfun(@(state) ~isempty(state), states);
    states = states(ind);
    
    time = cellfun(@(state) state.time, states);
    
    % Store results for this simulation
    all_states{sim_idx} = states;
    all_time{sim_idx} = time;
    
    % Evaluate lithium plating variables
    lp = 'LithiumPlating';
    
    varsToEval = { ...
        {'LithiumPlating', 'surfaceCoverage'}, ...
        {'LithiumPlating', 'platingFlux'}, ...
        {'LithiumPlating', 'chemicalFlux'}, ...
        {'LithiumPlating', 'etaPlating'}, ...
        {'LithiumPlating', 'etaChemical'} ...
                 };
    for k = 1:numel(states)
        for var = 1:numel(varsToEval)
            states{k} = model.evalVarName(states{k}, varsToEval{var});
        end
    end
    
    % Extract variables for plotting
    results{sim_idx}.time = time;
    results{sim_idx}.surfaceCoverage = cellfun(@(s) s.(lp).surfaceCoverage, states);
    results{sim_idx}.platingFlux = cellfun(@(s) s.(lp).platingFlux, states);
    results{sim_idx}.chemicalFlux = cellfun(@(s) s.(lp).chemicalFlux, states);
    results{sim_idx}.platedConcentration = cellfun(@(s) s.(lp).platedConcentration, states);
    results{sim_idx}.etaPlating = cellfun(@(s) s.(lp).etaPlating, states);
    results{sim_idx}.potentialWithoutActivity = cellfun(@(s) s.(lp).phiElectrode - s.(lp).phiElectrolyte, states);
    results{sim_idx}.kInter = kInter_values(sim_idx);
    
    fprintf('Simulation %d/%d completed.\n', sim_idx, n_simulations);
end

%% Plotting results with annotations

% Define colors for different simulations
colors = {'b', 'r', 'g'};
line_styles = {'-', '--', ':'};

% Plot 1: Surface Coverage
figure('Position', [100, 100, 800, 600]);
hold on;
for sim_idx = 1:n_simulations
    plot(results{sim_idx}.time/hour, results{sim_idx}.surfaceCoverage, ...
         'Color', colors{sim_idx}, 'LineStyle', line_styles{sim_idx}, 'LineWidth', 2, ...
         'DisplayName', sprintf('kInter = %.0e', results{sim_idx}.kInter));
end
xlabel('Time [hour]');
ylabel('Surface Coverage');
title('Surface Coverage of Plated Lithium for Different kInter Values');
legend('Location', 'best');
grid on;
hold off;

% Plot 2: Plating Flux
figure('Position', [200, 100, 800, 600]);
hold on;
for sim_idx = 1:n_simulations
    plot(results{sim_idx}.time/hour, results{sim_idx}.platingFlux, ...
         'Color', colors{sim_idx}, 'LineStyle', line_styles{sim_idx}, 'LineWidth', 2, ...
         'DisplayName', sprintf('kInter = %.0e', results{sim_idx}.kInter));
end
xlabel('Time [hour]');
ylabel('Plating Flux [mol/m²/s]');
title('Lithium Plating Flux for Different kInter Values');
legend('Location', 'best');
grid on;
hold off;

% Plot 3: Plated Concentration
figure('Position', [300, 100, 800, 600]);
hold on;
for sim_idx = 1:n_simulations
    plot(results{sim_idx}.time/hour, results{sim_idx}.platedConcentration, ...
         'Color', colors{sim_idx}, 'LineStyle', line_styles{sim_idx}, 'LineWidth', 2, ...
         'DisplayName', sprintf('kInter = %.0e', results{sim_idx}.kInter));
end
xlabel('Time [hour]');
ylabel('Plated Concentration [mol/m²]');
title('Accumulated Plated Lithium for Different kInter Values');
legend('Location', 'best');
grid on;
hold off;

% Plot 4: Eta Plating
figure('Position', [400, 100, 800, 600]);
hold on;
for sim_idx = 1:n_simulations
    plot(results{sim_idx}.time/hour, results{sim_idx}.etaPlating, ...
         'Color', colors{sim_idx}, 'LineStyle', line_styles{sim_idx}, 'LineWidth', 2, ...
         'DisplayName', sprintf('kInter = %.0e', results{sim_idx}.kInter));
end
xlabel('Time [hour]');
ylabel('Eta Plating [V]');
title('Plating Overpotential for Different kInter Values');
legend('Location', 'best');
grid on;
hold off;

% Plot 5: Potential without Activity
figure('Position', [500, 100, 800, 600]);
hold on;
for sim_idx = 1:n_simulations
    plot(results{sim_idx}.time/hour, results{sim_idx}.potentialWithoutActivity, ...
         'Color', colors{sim_idx}, 'LineStyle', line_styles{sim_idx}, 'LineWidth', 2, ...
         'DisplayName', sprintf('kInter = %.0e', results{sim_idx}.kInter));
end
xlabel('Time [hour]');
ylabel('Potential [V]');
title('Potential Without Activity for Different kInter Values');
legend('Location', 'best');
grid on;
hold off;
