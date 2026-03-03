%% run stand-alone active material model

% clear the workspace and close open figures
clear
close all

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa

%% Setup the properties of Li-ion battery materials and cell design
jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));

% We define some shorthand names for simplicity.
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

jsonstruct.(ne).(co).(am).isRootSimulationModel = true;

inputparams = BatteryInputParams(jsonstruct);

inputparams = inputparams.(ne).(co).(am);

model = ActiveMaterial(inputparams);
model = model.setupForSimulation();

%% Setup initial state

% shortcuts

sd  = 'SolidDiffusion';
itf = 'Interface';

cElectrolyte   = 5e-1*mol/litre;
phiElectrolyte = 0;
T              = 298;

cElectrodeInit   = (model.(itf).guestStoichiometry100)*(model.(itf).saturationConcentration);

% set primary variables
N = model.(sd).N;
initState.(sd).c        = cElectrodeInit*ones(N, 1);
initState.(sd).cSurface = cElectrodeInit;

% set static variable fields
initState.T = T;
initState.(itf).cElectrolyte   = cElectrolyte;
initState.(itf).phiElectrolyte = phiElectrolyte;

initState = model.evalVarName(initState, {itf, 'OCP'});

OCP = initState.(itf).OCP;
initState.E = OCP + phiElectrolyte;

%% setup schedule

% Reference rate which roughly corresponds to 1 hour for the data of this example
Iref = 5e-12;

Imax = 2e1*Iref;

total = 1*hour*(Iref/Imax);
n     = 100;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

tup = 1*second*(Iref/Imax); % rampup value for the current function, see rampupSwitchControl
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

figure
plot(time/hour, E);
xlabel('time [hour]');
ylabel('Potential [mol/L]');

cmin = cellfun(@(state) min(state.(sd).c), states);
cmax = cellfun(@(state) max(state.(sd).c), states);

for istate = 1 : numel(states)
    states{istate} = model.evalVarName(states{istate}, {sd, 'cAverage'});
end

caver = cellfun(@(state) max(state.(sd).cAverage), states);

figure
hold on
plot(time/hour, cmin /(mol/litre), 'displayname', 'cmin');
plot(time/hour, cmax /(mol/litre), 'displayname', 'cmax');
plot(time/hour, caver/(mol/litre), 'displayname', 'total concentration');
title('Concentration in particle / mol/L')
legend show

c = states{end}.(sd).c;
r = linspace(0, model.(sd).particleRadius, model.(sd).N);

figure
plot(r, c/(mol/litre));
xlabel('radius / m')
ylabel('concentration / mol/L')
title('Particle concentration profile (last time step)')



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
