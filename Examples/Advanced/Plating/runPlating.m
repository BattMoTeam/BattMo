%% run stand-alone active material model with lithium plating

clear
close all

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

jsonstruct.(ne).(co).(am).useLithiumPlating = false;

% Flag pour modèle stand-alone
jsonstruct.(ne).(co).(am).isRootSimulationModel = true;

jsonstruct_lithium_plating = parseBattmoJson(fullfile('Examples', 'Advanced', 'Plating', 'lithium_plating.json'));

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


F = model.(itf).constants.F;
R = model.(itf).constants.R;

if model.useLithiumPlating
    nPl0                 = model.LithiumPlating.nPl0;
    r                    = model.LithiumPlating.particleRadius;
    vf                   = model.LithiumPlating.volumeFraction;
    platedConcentration0 = nPl0 * vf / ((4/3)*pi*r^3);
    
    %initialisation so that overpotential are =0 at the beginning
    platedConcentrationInit = platedConcentration0/(exp((F*OCP)/(R*T)) - 1)^(1/4);

    model.(lp).platedConcentrationRef = platedConcentrationInit;

    initState.(lp).platedConcentrationNorm = platedConcentrationInit / model.(lp).platedConcentrationRef;
    initState.(lp).phiSolid            = initState.E;
    initState.(lp).phiElectrolyte      = phiElectrolyte;
    initState.(lp).cElectrolyte        = cElectrolyte;
    initState.(lp).nSEI                = 0;
end

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
    cmax = (model.(itf).guestStoichiometry100)*(model.(itf).saturationConcentration);
    control.stopFunction = @(model, state, state0_inner) (state.(sd).cSurface >= cmax);
  case 'discharge'
    srcfunc = @(time) rampupControl(time, tup, Imax); %0 pour tourner à vide
    cmin = (model.(itf).guestStoichiometry0)*(model.(itf).saturationConcentration);
    control.stopFunction = @(model, state, state0_inner) (state.(sd).cSurface <= cmin);
  otherwise
    error('scenario not recognized');
end

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
nls.errorOnFailure = false;
model.nonlinearTolerance = 1e-2;

%% Run simulation
% model.G = cartGrid([1, 1]);  % grille fictive 1x1 si pas de spatialisation
% model.G = computeGeometry(model.G); !!!

model.verbose = true;
[~, states, report] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

%% plotting

ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

time     = cellfun(@(state) state.time, states);
cSurface = cellfun(@(state) state.(sd).cSurface, states);
E        = cellfun(@(state) state.E, states);

figure
plot(time/hour, cSurface/(1/litre));
xlabel('time [hour]');
ylabel('Surface concentration [mol/L]');
title('Surface concentration');

figure
plot(time/hour, E);
xlabel('time [hour]');
ylabel('Potential [mol/L]');
title('Potential difference');

cmin = cellfun(@(state) min(state.(sd).c), states);
cmax = cellfun(@(state) max(state.(sd).c), states);

for istate = 1 : numel(states)
    states{istate} = model.evalVarName(states{istate}, {sd, 'cAverage'});
end

caver = cellfun(@(state) max(state.(sd).cAverage), states);

% figure
% hold on
% plot(time/hour, cmin /(mol/litre), 'displayname', 'cmin');
% plot(time/hour, cmax /(mol/litre), 'displayname', 'cmax');
% plot(time/hour, caver/(mol/litre), 'displayname', 'total concentration');
% title('Concentration in particle / mol/L')
% xlabel('time [hour]');
% ylabel('Concentration [mol/L]');
% legend show
% 
% c = states{end}.(sd).c;
% r = linspace(0, model.(sd).particleRadius, model.(sd).N);
% 
% figure
% plot(r, c/(mol/litre));
% xlabel('radius / m')
% ylabel('concentration / mol/L')
% title('Particle concentration profile (last time step)')

%% Lithium plating plotting

lp = 'LithiumPlating';

varsToEval = {{'Interface'     , 'eta'}         , ...
              {'LithiumPlating', 'etaPlating'}  , ...
              {'LithiumPlating', 'etaChemical'} , ...
              {'Interface'     , 'intercalationFlux'}           , ...
              {'LithiumPlating', 'platingFlux'} , ...
              {'LithiumPlating', 'chemicalFlux'}, ...
              {'LithiumPlating', 'surfaceCoverage'}};
for k = 1:numel(states)
    for var = 1:numel(varsToEval)
        states{k} = model.evalVarName(states{k}, varsToEval{var});
    end
end

varnames = {'eta', ...             
            'etaPlating', ...      
            'etaChemical', ...     
            'platingFlux', ...     
            'chemicalFlux', ...    
            'intercalationFlux', ...               
            'surfaceCoverage', ...
            'platedConcentrationNorm'};

vars = {};

vars{end + 1} = cellfun(@(s) s.(itf).eta, states);
vars{end + 1} = cellfun(@(s) s.(lp).etaPlating, states);
vars{end + 1} = cellfun(@(s) s.(lp).etaChemical, states);
vars{end + 1} = cellfun(@(s) s.(lp).platingFlux, states);
vars{end + 1} = cellfun(@(s) s.(lp).chemicalFlux, states);
vars{end + 1} = cellfun(@(s) s.(itf).intercalationFlux, states);
vars{end + 1} = cellfun(@(s) s.(lp).surfaceCoverage, states);
vars{end + 1} = cellfun(@(s) s.(lp).platedConcentrationNorm, states);

for ivar = 1 : numel(varnames)
    figure
    plot(time, vars{ivar});
    title(varnames{ivar});
end

