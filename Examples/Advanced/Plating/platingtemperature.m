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

%% Simulation setup for multiple temperatures

T_list = [150, 250, 350];  % différentes températures
styles = {'-r', '-g', '-b'};  % couleurs et pointillés
results = struct();

Iref = 3e-13;
Imax = Iref;
total = 1e-2*hour*(Iref/Imax);
n = 5000;
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

% Variables à évaluer
varsToEval = {{'Interface','eta'}, ...
              {'LithiumPlating','etaPlating'}, ...
              {'LithiumPlating','etaChemical'}, ...
              {'LithiumPlating','platingFlux'}, ...
              {'LithiumPlating','chemicalFlux'}, ...
              {'Interface','intercalationFlux'}, ...
              {'LithiumPlating','surfaceCoverage'}, ...
              {'LithiumPlating','platedConcentration'}};

for iT = 1:length(T_list)
    T = T_list(iT);
    
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
    time = cellfun(@(state) state.time, states);
    cSurface = cellfun(@(state) state.(sd).cSurface, states);
    E = cellfun(@(state) state.E, states);
    
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
    results(iT).caver = caver;
    results(iT).lp_vars = lp_vars;
end

%% Plotting all variables

% varnames = {'eta', 'etaPlating', 'etaChemical', ...
%             'platingFlux', 'chemicalFlux', 'intercalationFlux', ...
%             'surfaceCoverage', 'platedConcentration'};
% 
% for ivar = 1:length(varnames)
%     figure; hold on;
%     for iT = 1:length(T_list)
%         plot(results(iT).time/hour, results(iT).lp_vars{ivar}, ...
%              styles{iT}, 'DisplayName', sprintf('T = %d K', T_list(iT)));
%     end
%     xlabel('Time [hour]');
%     ylabel(varnames{ivar});
%     title(varnames{ivar});
%     legend show;
% end
% 
% figure; hold on;
% for iT = 1:length(T_list)
%     plot(results(iT).time/hour, results(iT).cSurface/(1/litre), styles{iT}, ...
%          'DisplayName', sprintf('T = %d K', T_list(iT)));
% end
% xlabel('Time [hour]');
% ylabel('Surface concentration [mol/L]');
% title('Surface concentration');
% legend show;
% 
% figure; hold on;
% for iT = 1:length(T_list)
%     plot(results(iT).time/hour, results(iT).E, styles{iT}, ...
%          'DisplayName', sprintf('T = %d K', T_list(iT)));
% end
% xlabel('Time [hour]');
% ylabel('Potential [V]');
% title('Potential difference');
% legend show;

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
h_c = gobjects(1, length(T_list));  % solid lines (concentration)
h_s = gobjects(1, length(T_list));  % dashed lines (surface coverage)

yyaxis left;
ylabel('Surface concentration [mol/L]');
ax = gca;
ax.YColor = 'k';  % black axis

for iT = 1:length(T_list)
    h_c(iT) = plot(results(iT).time/hour, ...
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

for iT = 1:length(T_list)
    surfaceCoverage = results(iT).lp_vars{7};  % surface coverage
    h_s(iT) = plot(results(iT).time/hour, ...
        surfaceCoverage, '--', ...
        'Color', customColors(iT,:), ...
        'LineWidth', 1.0, ...
        'LineStyle', '--');
end

xlabel('Time [hour]');
title('Surface concentration and surface coverage vs time');
grid on;

% Build custom legend
% One entry per color (temperature)
temp_legend = gobjects(1, length(T_list));
for iT = 1:length(T_list)
    temp_legend(iT) = plot(nan, nan, '-', 'Color', customColors(iT,:), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('T = %d K', T_list(iT)));
end

% One entry per line style (solid = concentration, dashed = coverage)
type_legend(1) = plot(nan, nan, '-', 'Color', 'k', ...
    'LineWidth', 1.5, 'DisplayName', 'Surface concentration');
type_legend(2) = plot(nan, nan, '--', 'Color', 'k', ...
    'LineWidth', 1.0, 'DisplayName', 'Surface coverage');

% Show combined legend
legend([temp_legend, type_legend], 'Location', 'best');

% Save figure as high-resolution PNG
saveas(gcf, 'surface_concentration_coverage_temperature.png');