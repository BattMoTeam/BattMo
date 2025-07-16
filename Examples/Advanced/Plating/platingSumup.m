%% Run stand-alone active material model with lithium plating
% This example shows how to simulate a single particle of a silicon graphite
% electrode, taking into account the plating phenomenon

clear all
close all

%% Setup the properties of Li-ion battery materials and cell design

jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));

inputparams = BatteryInputParams(jsonstruct);

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

% Flag for stand-alone model
jsonstruct.(ne).(co).(am).isRootSimulationModel = true;

jsonstruct_lithium_plating = parseBattmoJson(fullfile('Examples', 'Advanced', 'Plating', 'lithium_plating.json'));

jsonstruct.(ne).(co).(am).LithiumPlating = jsonstruct_lithium_plating.LithiumPlating;

scenario = 'charge';

%% following is not used at particle level

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

%OCP is computed via a function described in the article (S-5)
inputparams.Interface.openCircuitPotential.functionname = 'computeOCP_Graphite_Latz';

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
    % cElectrodeInit = (model.(itf).guestStoichiometry0)*(model.(itf).saturationConcentration);
    cElectrodeInit = 29.5*mol/litre;
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


F = model.(itf).constants.F;
R = model.(itf).constants.R;

if model.useLithiumPlating
    nPl0                 = model.LithiumPlating.nPl0;
    r                    = model.LithiumPlating.particleRadius;
    vf                   = model.LithiumPlating.volumeFraction;
    platedConcentration0 = nPl0 * vf / ((4/3)*pi*r^3);
    
    %initialisation so that overpotential are = 0 and reaction at equilibrium
    %at the beginning
    platedConcentrationInit = platedConcentration0/(exp((F*OCP)/(R*T)) - 1)^(1/4);

    model.(lp).platedConcentrationRef = platedConcentrationInit;

    initState.(lp).platedConcentrationNorm = platedConcentrationInit / model.(lp).platedConcentrationRef;
    initState.(lp).phiSolid            = initState.E;
    initState.(lp).phiElectrolyte      = phiElectrolyte;
    initState.(lp).cElectrolyte        = cElectrolyte;
    initState.(lp).nSEI                = 0;
end

%% setup schedule
%This part is essential to see the lithium plating effect. Increasing Iref
%strengthen the effect.

Iref = 7e-13;
Imax = Iref;
total = 1e-1*hour*(Iref/Imax); %total time of the charge
n     = 500;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

tup = 1*second*(Iref/Imax);

switch scenario
  case 'charge'
    srcfunc = @(time) rampupControl(time, tup, -Imax);
    cmax = (model.(itf).guestStoichiometry100)*(model.(itf).saturationConcentration);
    % control.stopFunction = @(model, state, state0_inner) (state.(sd).cSurface >= cmax);
  case 'discharge'
    srcfunc = @(time) rampupControl(time, tup, Imax);
    cmin = (model.(itf).guestStoichiometry0)*(model.(itf).saturationConcentration);
    control.stopFunction = @(model, state, state0_inner) (state.(sd).cSurface <= cmin);
  otherwise
    error('scenario not recognized');
end

% uncomment to make the simulation run without current (debugging)
% srcfunc = @(time) 0;

control.src = srcfunc;

schedule = struct('control', control, 'step', step);

scalingparams = struct('I'                  , Imax                              , ...
                       'elyteConcentration' , initState.(itf).cElectrolyte);

if model.useLithiumPlating
    scalingparams.platedConcentration = platedConcentrationInit;
end


model = model.setupScalings(scalingparams);

%% setup non-linear solver

nls = NonLinearSolver();
nls.errorOnFailure  = false;
nls.maxTimestepCuts = 20;

model.nonlinearTolerance = 1e-6;

%% Run simulation by charging the particle

model.verbose = true;
[~, states, report] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

%% Then discharge it
%Comment this part to only see the charge

%We use the last state of the previous simulation to initialise the new one
initstate = states{end};
jsonstruct.(ctrl).DRate = -1;
jsonstruct.Control.controlPolicy = 'CCDischarge';

% And we save the charging states for later

chargeStates = states;

% Setup schedule

% Set up a StopFunction if necessary
% Control.stopFunction = @(model, state, state0_inner) (1 == 0);

srcfunc = @(time) rampupControl(time, tup, Imax);

control.src = srcfunc;

schedule = struct('control', control, 'step', step);

% Run simulation

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

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

cmin = cellfun(@(state) min(state.(sd).c), states);
cmax = cellfun(@(state) max(state.(sd).c), states);

for istate = 1 : numel(states)
    states{istate} = model.evalVarName(states{istate}, {sd, 'cAverage'});
end


%% Simple plots 
% Necessary variables
for k = 1:numel(states)
    states{k} = model.evalVarName(states{k}, {'Interface', 'intercalationFlux'});
    states{k} = model.evalVarName(states{k}, {'LithiumPlating', 'platingFlux'});
    states{k} = model.evalVarName(states{k}, {'LithiumPlating', 'surfaceCoverage'});
end

% flux extraction
platingFlux = cellfun(@(s) s.LithiumPlating.platingFlux .* s.LithiumPlating.surfaceCoverage, states);
intercalationFlux = cellfun(@(s) s.Interface.intercalationFlux .* (1 - s.LithiumPlating.surfaceCoverage), states);

% Tracé des deux courbes
figure
plot(time, platingFlux, '-', 'DisplayName', 'Plating Flux'); hold on
plot(time, intercalationFlux, '-', 'DisplayName', 'Intercalation Flux');
xlabel('Time [second]');
ylabel('Flux [mol/m²/s]');
title('Plating vs Intercalation Flux');
legend show
grid on