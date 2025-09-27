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
jsonstruct.(pe).(co).(am).diffusionModelType = 'full';

jsonstruct.(ne).(co).(am).useLithiumPlating = true;

%%
% OCP is computed via a function described in the article (S-5)
jsonstruct.(ne).(co).(am).(itf).openCircuitPotential.functionName = 'computeOCP_Graphite_Latz';
jsonstruct.(ne).(co).(am).(itf).includeEntropyChange = false;


%%
% Flag for stand-alone model
%

jsonstruct.(ne).(co).(am).isRootSimulationModel = true;

jsonstruct_lithium_plating = parseBattmoJson(fullfile('Examples', 'Advanced', 'Plating', 'lithium_plating.json'));

jsonstruct.(ne).(co).(am).LithiumPlating = jsonstruct_lithium_plating.LithiumPlating;

jsonstruct.NegativeElectrode.Coating.ActiveMaterial.LithiumPlating.kPl = 1e2*jsonstruct.NegativeElectrode.Coating.ActiveMaterial.LithiumPlating.kInter;


%%
% Setup InputParams

inputparams = BatteryInputParams(jsonstruct);
inputparams = inputparams.(ne).(co).(am);


%% Setup the model

model = ActiveMaterial(inputparams);

%%
% We equip the model for simulation
%

model = model.setupForSimulation();
model.verbose = true;



%% setup schedule
%This part is essential to see the lithium plating effect. Increasing Iref
%strengthen the effect.

Iref = 7e-13;
Imax = Iref;
total = 3e-2*hour*(Iref/Imax); %total time of the charge
n     = 500;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

tup = 1*second*(Iref/Imax);

srcfunc = @(time) rampupControl(time, tup, -Imax);
control.src = srcfunc;

schedule = struct('control', control, 'step', step);


%% Setup initial state

sd  = 'SolidDiffusion';
itf = 'Interface';

cElectrolyte   = 5e-1*mol/litre;
phiElectrolyte = 0;
T              = 298;
cElectrodeInit = 29.9*mol/litre;

[model, initstate] = setupPlatingInitialState(model, T, cElectrolyte, phiElectrolyte, cElectrodeInit, Imax);


%% setup non-linear solver

nls = NonLinearSolver();
nls.errorOnFailure  = false;
nls.maxTimestepCuts = 20;

model.nonlinearTolerance = 1e-6;

%% Run simulation by charging the particle

inputSim = struct('model'          , model    , ...
                  'schedule'       , schedule , ...
                  'initstate'      , initstate, ...
                  'NonLinearSolver', nls);
simsetup = SimulationSetup(inputSim);

states = simsetup.run();

chargeStates = states; % for later


initstate = states{end};

srcfunc = @(time) rampupControl(time, tup, Imax);
control.src = srcfunc;
schedule = struct('control', control, 'step', step);


simsetup.initstate = initstate;
simsetup.schedule = schedule;

% Run simulation

states = simsetup.run();

% And concatenate the states

dischargeStates = states;
states = vertcat(chargeStates, dischargeStates);

%% Plot

ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

time     = cellfun(@(state) state.time, states);
cSurface = cellfun(@(state) state.(sd).cSurface, states);
E        = cellfun(@(state) state.E, states);

figure
plot(time, cSurface/(1/litre));
xlabel('time [second]');
ylabel('Surface concentration [mol/L]');
title('Surface concentration');

figure
plot(time, E);
xlabel('time [second]');
ylabel('Potential [mol/L]');
title('Potential difference');


%% plots
%
% Update the flux variables
for istate = 1:numel(states)
    states{istate} = model.evalVarName(states{istate}, {sd, 'cAverage'});
    states{istate} = model.evalVarName(states{istate}, {'Interface', 'intercalationFlux'});
    states{istate} = model.evalVarName(states{istate}, {'LithiumPlating', 'platingFlux'});
    states{istate} = model.evalVarName(states{istate}, {'LithiumPlating', 'surfaceCoverage'});
end

vsa = model.(lp).volumetricSurfaceArea;

%%
% Retrieve the fluxes
platingFlux       = cellfun(@(s) s.LithiumPlating.platingFlux .* s.LithiumPlating.surfaceCoverage * vsa, states);
intercalationFlux = cellfun(@(s) s.Interface.intercalationFlux .* (1 - s.LithiumPlating.surfaceCoverage ) * vsa, states);

%%
% Plot of the two fluxes
figure
plot(time, platingFlux, '-', 'DisplayName', 'Plating Flux'); hold on
plot(time, intercalationFlux, '-', 'DisplayName', 'Intercalation Flux');
xlabel('Time [second]');
ylabel('Volumetric Flux [mol/m³/s]');
title('Plating vs Intercalation Flux');
legend show
grid on


%% Plot flux as percentage of Imax

F = model.(itf).constants.F; % Faraday constant (C/mol)
n = 1; % numner of exchanged electrons (1 for Li^+)

vsa = model.(lp).volumetricSurfaceArea;

platingFlux       = cellfun(@(s) s.LithiumPlating.platingFlux.*s.LithiumPlating.surfaceCoverage*vsa, states);
intercalationFlux = cellfun(@(s) s.Interface.intercalationFlux.*(1 - s.LithiumPlating.surfaceCoverage)*vsa, states);

platingCurrent       = platingFlux*n*F;
intercalationCurrent = intercalationFlux*n*F;

totalCurrent = abs(platingCurrent + intercalationCurrent);

platingPct       = 100*platingCurrent./totalCurrent;
intercalationPct = 100*intercalationCurrent./totalCurrent;

%%
% Plot
figure
plot(time, platingPct, '-', 'LineWidth', 1.5, 'DisplayName', 'Plating'); hold on
plot(time, intercalationPct, '-', 'LineWidth', 1.5, 'DisplayName', 'Intercalation');
xlabel('Time [s]');
ylabel('Current [% of I_{max}]');
title('Plating vs Intercalation Current as % of I_{max}');
legend show
grid on

