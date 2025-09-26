%% Run stand-alone active material model with lithium plating
% This example shows how to simulate a single particle of a silicon graphite
% electrode, taking into account the plating phenomenon

clear all
close all

%% Setup the properties of Li-ion battery materials and cell design

filename = fullfile('ParameterData'        , ...
                    'BatteryCellParameters', ...
                    'LithiumIonBatteryCell', ...
                    'lithium_ion_battery_nmc_graphite.json');

jsonstruct = parseBattmoJson(filename);

%%
% We define some shortcuts

ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
lp      = 'LithiumPlating';
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
jsonstruct.(ne).(co).(am).useLithiumPlating  = true;

% OCP is computed via a function described in the article (S-5)
jsonstruct.(ne).(co).(am).(itf).openCircuitPotential.functionName = 'computeOCP_Graphite_Latz';

%%
% Add lithium plating parameters
%

jsonstruct_lithium_plating = parseBattmoJson(fullfile('Examples', 'Advanced', 'Plating', 'lithium_plating.json'));

jsonstruct.(ne).(co).(am).LithiumPlating = jsonstruct_lithium_plating.LithiumPlating;


%% Setup the model

%%
% Flag for stand-alone model
%

jsonstruct.(ne).(co).(am).isRootSimulationModel = true;

inputparams = BatteryInputParams(jsonstruct);
inputparams = inputparams.(ne).(co).(am);

model = ActiveMaterial(inputparams);

%%
% We equip the model for simulation
%

model = model.setupForSimulation();
model.verbose = true;

%% Setup the schedule

Iref  = 3e-13;
Imax  = Iref;
total = 1e-2*hour*(Iref/Imax);
n     = 500;
n     = 100;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));
tup   = 1*second*(Iref/Imax);

srcfunc = @(time) rampupControl(time, tup, -Imax);

control.src = srcfunc;
schedule = struct('control', control, 'step', step);

nls                      = NonLinearSolver();
nls.errorOnFailure       = false;
nls.maxTimestepCuts      = 20;
model.nonlinearTolerance = 1e-6;

inputSimSetup = struct('model'          , model, ...
                       'NonLinearSolver', nls  , ...
                       'schedule'       , schedule);


simsetup = SimulationSetup(inputSimSetup);


%%
% Variables to evaluate
%

varsToEval = {{'Interface'     ,'eta'}              , ...
              {'LithiumPlating','etaPlating'}       , ...
              {'LithiumPlating','etaChemical'}      , ...
              {'LithiumPlating','platingFlux'}      , ...
              {'LithiumPlating','chemicalFlux'}     , ...
              {'Interface'     ,'intercalationFlux'}, ...
              {'LithiumPlating','surfaceCoverage'}  , ...
              {'LithiumPlating','platedConcentration'}};

%%
% Initial values, common to all simulations
%

cElectrolyte   = 5e-1*mol/litre;
phiElectrolyte = 0;
cElectrodeInit = 30*mol/litre;


%% Run first simulation, particle charge
%

T_list = [150, 250, 350];  % different temperatures
T_list = [350];  % different temperatures
styles = {'-r', '-g', '-b'};  % couleurs et pointillés
results = struct();

for iT = 1:length(T_list)
    
    T = T_list(iT);

    model = simsetup.model;
    [model, initstate] = setupPlatingInitialState(model, ...
                                                  T             , ...
                                                  cElectrolyte  , ...
                                                  phiElectrolyte, ...
                                                  cElectrodeInit, ...
                                                  Imax);
    simsetup.model     = model;
    simsetup.initstate = initstate;
    

    %% Run simulation
    %

    states = simsetup.run();

    return
    
    states = states(cellfun(@(s) ~isempty(s), states));
    
    time     = cellfun(@(state) state.time, states);
    cSurface = cellfun(@(state) state.(sd).cSurface, states);
    E        = cellfun(@(state) state.E, states);
    
    for istate = 1:numel(states)
        states{istate} = simsetup.model.evalVarName(states{istate}, {sd, 'cAverage'});
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

    results(iT).T        = T;
    results(iT).time     = time;
    results(iT).cSurface = cSurface;
    results(iT).E        = E;
    results(iT).caver    = caver;
    results(iT).lp_vars  = lp_vars;
    
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
