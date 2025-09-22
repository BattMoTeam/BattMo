%% run stand-alone active material model with lithium plating

clear all
close all

%% Setup the properties of Li-ion battery materials and cell design

jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
inputparams = BatteryInputParams(jsonstruct);

use_cccv = true;
if use_cccv
    inputparams.SOC = 0;
    cccvstruct = struct('controlPolicy', 'CCCV', ...
                        'initialControl', 'charging', ...
                        'numberOfCycles', 2, ...
                        'CRate', 1.5, ...
                        'DRate', 1, ...
                        'lowerCutoffVoltage', 3, ...
                        'upperCutoffVoltage', 4, ...
                        'dIdtLimit', 1e-2, ...
                        'dEdtLimit', 1e-4);
    cccvinputparams = CcCvControlModelInputParams(cccvstruct);
    inputparams.Control = cccvinputparams;
end

ne = 'NegativeElectrode';
pe = 'PositiveElectrode';
lp = 'LithiumPlating';
elyte = 'Electrolyte';
thermal = 'ThermalModel';
co = 'Coating';
am = 'ActiveMaterial';
itf = 'Interface';
sd = 'SolidDiffusion';
ctrl = 'Control';
cc = 'CurrentCollector';

jsonstruct.use_thermal = false;
jsonstruct.include_current_collectors = false;

jsonstruct.(ne).(co).(am).diffusionModelType = 'full';
jsonstruct.(pe).(co).(am).diffusionModelType = 'full';
jsonstruct.(ne).(co).(am).useLithiumPlating = true;
jsonstruct.(ne).(co).(am).isRootSimulationModel = true;

jsonstruct_lithium_plating = parseBattmoJson(fullfile('Examples', 'Advanced', 'Plating', 'lithium_plating.json'));
jsonstruct.(ne).(co).(am).LithiumPlating = jsonstruct_lithium_plating.LithiumPlating;

scenario = 'charge';
switch scenario
    case 'charge'
        jsonstruct.Control.controlPolicy = 'CCCharge';
    case 'discharge'
        jsonstruct.Control.controlPolicy = 'CCDischarge';
    otherwise
        error('scenario not recognized');
end

inputparams = BatteryInputParams(jsonstruct);
inputparams = inputparams.(ne).(co).(am);
inputparams.Interface.openCircuitPotential.functionname = 'computeOCP_Graphite_Latz';
model = ActiveMaterial(inputparams);
model = model.setupForSimulation();

%% Simulation setup for multiple intensities

I_list = [1e-13, 3e-13, 5e-13];  % différentes intensities
styles = {'-r', '-g', '-b'};  % couleurs et pointillés
results = struct();

T = 300;

% Variables à évaluer
varsToEval = {{'Interface','eta'}, ...
              {'LithiumPlating','etaPlating'}, ...
              {'LithiumPlating','etaChemical'}, ...
              {'LithiumPlating','platingFlux'}, ...
              {'LithiumPlating','chemicalFlux'}, ...
              {'Interface','intercalationFlux'}, ...
              {'LithiumPlating','surfaceCoverage'}, ...
              {'LithiumPlating','platedConcentration'}};

for iT = 1:length(I_list)
    Imax = I_list(iT);
    Iref = Imax;
    total = 2e-2*hour*(Iref/Imax);
    n = 500;
    dt = total/n;
    step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));
    tup = 1*second*(Iref/Imax);

    switch scenario
        case 'charge'
            srcfunc = @(time) rampupControl(time, tup, -Imax);
        case 'discharge'
            srcfunc = @(time) rampupControl(time, tup, Imax);
        otherwise
            error('scenario not recognized');
    end

    control.src = srcfunc;
    schedule = struct('control', control, 'step', step);

    nls = NonLinearSolver();
    nls.errorOnFailure  = false;
    nls.maxTimestepCuts = 20;
    model.nonlinearTolerance = 1e-6;

    %% Setup initial state
    cElectrolyte = 5e-1*mol/litre;
    phiElectrolyte = 0;

    switch scenario
        case 'charge'
            cElectrodeInit = 30*mol/litre;
        case 'discharge'
            cElectrodeInit = (model.(itf).guestStoichiometry100)*(model.(itf).saturationConcentration);
    end

    N = model.(sd).N;
    initState.(sd).c = cElectrodeInit*ones(N, 1);
    initState.(sd).cSurface = cElectrodeInit;
    initState.T = T;
    initState.(itf).cElectrolyte = cElectrolyte;
    initState.(itf).phiElectrolyte = phiElectrolyte;

    initState = model.evalVarName(initState, {itf, 'OCP'});
    OCP = initState.(itf).OCP;
    initState.E = OCP + phiElectrolyte;

    F = model.(itf).constants.F;
    R = model.(itf).constants.R;

    if model.useLithiumPlating
        nPl0 = model.LithiumPlating.nPl0;
        r = model.LithiumPlating.particleRadius;
        vf = model.LithiumPlating.volumeFraction;
        platedConcentration0 = nPl0 * vf / ((4/3)*pi*r^3);
        platedConcentrationInit = platedConcentration0/(exp((F*OCP)/(R*T)) - 1)^(1/4);

        model.(lp).platedConcentrationRef = platedConcentrationInit;
        initState.(lp).platedConcentrationNorm = platedConcentrationInit / model.(lp).platedConcentrationRef;
        initState.(lp).phiSolid = initState.E;
        initState.(lp).phiElectrolyte = phiElectrolyte;
        initState.(lp).cElectrolyte = cElectrolyte;
        initState.(lp).nSEI = 0;
    end

    scalingparams = struct('I', Imax, ...
                           'elyteConcentration', initState.(itf).cElectrolyte);
    if model.useLithiumPlating
        scalingparams.platedConcentration = platedConcentrationInit;
    end

    model = model.setupScalings(scalingparams);

    %% Run simulation
    model.verbose = false;
    [~, states, ~] = simulateScheduleAD(initState, model, schedule, ...
        'OutputMinisteps', true, 'NonLinearSolver', nls);
    
    states = states(cellfun(@(s) ~isempty(s), states));
    states = states(cellfun(@(s) ~isempty(s), states));
    states = cellfun(@(s) model.evalVarName(s, {sd, 'cAverage'}), states, 'UniformOutput', false);

    % Time
    time = cellfun(@(s) s.time, states);
    time = [0; time];

    % Applied current
    I_applied = arrayfun(srcfunc, time(2:end));

    % Capacity
    Q = -diff(time) .* I_applied;  % Q(t) = ∫ I(t) dt

    % Cumulated capacity [mAh]
    Q_mAh = cumsum(Q) / 3.6;

    % Surface concentration et potentiel
    cSurface = cellfun(@(s) s.(sd).cSurface, states);
    E = cellfun(@(s) s.E, states);
    for istate = 1:numel(states)
        states{istate} = model.evalVarName(states{istate}, {sd, 'cAverage'});
    end
    caver = cellfun(@(state) max(state.(sd).cAverage), states);

    for k = 1:numel(states)
        for var = 1:numel(varsToEval)
            states{k} = model.evalVarName(states{k}, varsToEval{var});
        end
    end

    lp_vars = {
        cellfun(@(s) s.(itf).eta, states);
        cellfun(@(s) s.(lp).etaPlating, states);
        cellfun(@(s) s.(lp).etaChemical, states);
        cellfun(@(s) s.(lp).platingFlux, states);
        cellfun(@(s) s.(lp).chemicalFlux, states);
        cellfun(@(s) s.(itf).intercalationFlux, states);
        cellfun(@(s) s.(lp).surfaceCoverage, states);
        cellfun(@(s) s.(lp).platedConcentration, states);
    };

    results(iT).T = T;
    results(iT).time = time;
    results(iT).cSurface = cSurface;
    results(iT).E = E;

    results(iT).Q_mAh = Q_mAh;

    results(iT).caver = caver;
    results(iT).lp_vars = lp_vars;
end

%% Plotting all variables in function of passed capacity (Q_mAh)

varnames = {'eta', 'etaPlating', 'etaChemical', ...
            'platingFlux', 'chemicalFlux', 'intercalationFlux', ...
            'surfaceCoverage', 'platedConcentration'};

for ivar = 1:length(varnames)
    figure; hold on;
    for iT = 1:length(I_list)
        Q_mAh = results(iT).Q_mAh;
        ydata = results(iT).lp_vars{ivar};
        plot(Q_mAh, ydata, styles{iT}, ...
             'DisplayName', sprintf('I = %.1e A', I_list(iT)));
    end
    xlabel('Capacity passed [mAh]');
    ylabel(varnames{ivar});
    title(sprintf('%s vs capacity passed', varnames{ivar}));
    legend show;
    grid on;
end

%Surface concentration
figure; hold on;
for iT = 1:length(I_list)
    Q_mAh = results(iT).Q_mAh;
    cSurf = results(iT).cSurface / (1/litre);
    plot(Q_mAh, cSurf, styles{iT}, ...
         'DisplayName', sprintf('I = %.1e A', I_list(iT)));
end
xlabel('Capacity passed [mAh]');
ylabel('Surface concentration [mol/L]');
title('Surface concentration vs capacity passed');
legend show;
grid on;

% Potential difference
figure; hold on;
for iT = 1:length(I_list)
    Q_mAh = results(iT).Q_mAh;
    E = results(iT).E;
    plot(Q_mAh, E, styles{iT}, ...
         'DisplayName', sprintf('I = %.1e A', I_list(iT)));
end
xlabel('Capacity passed [mAh]');
ylabel('Potential [V]');
title('Potential difference vs capacity passed');
legend show;
grid on;

%% Plot the final graph

%% Combined Plot: Surface concentration + Surface coverage

figure;
hold on;

% Define one color per temperature
customColors = [ ...
    0.0000, 0.4470, 0.7410;  % blue
    0.8500, 0.3250, 0.0980;  % orange
    0.4660, 0.6740, 0.1880   % green
];

% Prepare graphics handles
h_c = gobjects(1, length(I_list));  % solid lines (concentration)
h_s = gobjects(1, length(I_list));  % dashed lines (surface coverage)

yyaxis left;
ylabel('Surface concentration [mol/L]');
ax = gca;
ax.YColor = 'k';  % black axis

for iT = 1:length(I_list)
    h_c(iT) = plot(results(iT).Q_mAh, ...
        results(iT).cSurface/(1/litre), '-', ...
        'Color', customColors(iT,:), ...
        'LineWidth', 1.5);
end

% % Add critical concentration line
% yline(30.0113, ':k', 'c_{crit}', ...
%     'LineWidth', 1.2, ...
%     'FontSize', 10, ...
%     'LabelHorizontalAlignment', 'left', ...
%     'LabelVerticalAlignment', 'bottom');

yyaxis right;
ylabel('Surface coverage [–]');
ax = gca;
ax.YColor = 'k';  % black axis

for iT = 1:length(I_list)
    surfaceCoverage = results(iT).lp_vars{7};  % surface coverage
    h_c(iT) = plot(results(iT).Q_mAh, ...
        surfaceCoverage, '--', ...
        'Color', customColors(iT,:), ...
        'LineWidth', 1.0, ...
        'LineStyle', '--');
end

xlabel('Time [hour]');
title('Surface concentration and surface coverage vs time (different currents)');
grid on;

% Build custom legend
% One entry per intensity
current_legend = gobjects(1, length(I_list));
for iT = 1:length(I_list)
    current_legend(iT) = plot(nan, nan, '-', 'Color', customColors(iT,:), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('I = %.1e A', I_list(iT)));
end

% One entry per line style
type_legend(1) = plot(nan, nan, '-', 'Color', 'k', ...
    'LineWidth', 1.5, 'DisplayName', 'Surface concentration');
type_legend(2) = plot(nan, nan, '--', 'Color', 'k', ...
    'LineWidth', 1.0, 'DisplayName', 'Surface coverage');

legend([current_legend, type_legend], 'Location', 'best');

% Save figure as high-resolution PNG
saveas(gcf, 'surface_concentration_coverage_temperature.png');