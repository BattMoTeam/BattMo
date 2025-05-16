%% run stand-alone active material model

clear
close all

%% Import the required modules from MRST
mrstModule add ad-core mrst-gui mpfa

%% Setup the properties of Li-ion battery materials and cell design

jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));

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

jsonstruct.(ne).(co).(am).useLithiumPlating = false;

% Flag pour modèle stand-alone
jsonstruct.(ne).(co).(am).isRootSimulationModel = true;

% Setup InputParams
inputparams = BatteryInputParams(jsonstruct);
inputparams = inputparams.(ne).(co).(am);

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

cElectrodeInit = (model.(itf).guestStoichiometry100)*(model.(itf).saturationConcentration);
N = model.(sd).N;
initState.(sd).c        = cElectrodeInit*ones(N, 1);
initState.(sd).cSurface = cElectrodeInit;

initState.T = T;
initState.(itf).cElectrolyte   = cElectrolyte;
initState.(itf).phiElectrolyte = phiElectrolyte;

initState = model.evalVarName(initState, {itf, 'OCP'});
OCP = initState.(itf).OCP;
initState.E = OCP + phiElectrolyte;

if model.useLithiumPlating
    lp = 'LithiumPlating';
    initState.(lp).nPl            = 0;
    initState.(lp).phiSolid       = initState.E;
    initState.(lp).phiElectrolyte = phiElectrolyte;
    initState.(lp).cElectrolyte   = cElectrolyte;
end

%% setup schedule

Iref = 5e-12;
Imax = 5e1*Iref;
total = 1*hour*(Iref/Imax);
n     = 100;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

tup = 1*second*(Iref/Imax);
srcfunc = @(time) rampupControl(time, tup, Imax);

cmin = (model.(itf).guestStoichiometry0)*(model.(itf).saturationConcentration);
control.stopFunction = @(model, state, state0_inner) (state.(sd).cSurface <= cmin);
control.src = srcfunc;

schedule = struct('control', control, 'step', step);

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

if model.useLithiumPlating
    lp = 'LithiumPlating';

    varsToEval = { ...
        {'LithiumPlating', 'surfaceCoverage'}, ...
        {'LithiumPlating', 'platingFlux'}, ...
        {'LithiumPlating', 'chemicalFlux'}, ...
        {'LithiumPlating', 'etaPlating'}, ...
        {'LithiumPlating', 'etaChemical'} ...
    };
    for k = 1:numel(states) %!!!
        for iv = 1:numel(varsToEval)
            states{k} = model.evalVarName(states{k}, varsToEval{iv});
        end
    end

    surfaceCoverage = cellfun(@(s) s.(lp).surfaceCoverage, states);
    platingFlux     = cellfun(@(s) s.(lp).platingFlux, states);
    chemicalFlux    = cellfun(@(s) s.(lp).chemicalFlux, states);
    nPl             = cellfun(@(s) s.(lp).nPl, states);

    figure
    plot(time/hour, surfaceCoverage);
    xlabel('time [hour]');
    ylabel('Surface coverage');
    title('Surface Coverage of Plated Lithium');

    figure
    plot(time/hour, platingFlux*model.(itf).volumetricSurfaceArea);
    xlabel('time [hour]');
    ylabel('Volumetric Plating Flux [mol/m³/s]');
    title('Lithium Plating Flux');

    figure
    plot(time/hour, nPl);
    xlabel('time [hour]');
    ylabel('nPl [mol/m²]');
    title('Accumulated Plated Lithium');

end
