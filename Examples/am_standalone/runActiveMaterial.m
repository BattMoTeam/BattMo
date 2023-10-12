%% run stand-alone active material model

% clear the workspace and close open figures
clear all
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

paramobj = BatteryInputParams(jsonstruct);

paramobj = paramobj.(ne).(co).(am);

paramobj.standAlone = true;

paramobj = paramobj.validateInputParams();

model = ActiveMaterial(paramobj);

cgt = model.computationalGraph();

return
%% Setup initial state

% shortcuts

sd  = 'SolidDiffusion';
itf = 'Interface';

cElectrolyte     = 5e-1*mol/litre;
phiElectrolyte   = 0;
T                = 298;

cElectrodeInit   = (model.(itf).theta100)*(model.(itf).cmax);

% set primary variables
N = model.(sd).N;
initState.(sd).c        = cElectrodeInit*ones(N, 1);
initState.(sd).cSurface = cElectrodeInit;

% set static variable fields
initState.T = T;
initState.(itf).cElectrolyte   = cElectrolyte;
initState.(itf).phiElectrolyte = phiElectrolyte;

initState = model.updateConcentrations(initState);
initState = model.dispatchTemperature(initState);
initState.(itf) = model.(itf).updateOCP(initState.(itf));

OCP = initState.(itf).OCP;
initState.phi = OCP + phiElectrolyte;

%% setup schedule

controlsrc = 1e3;

total = (10*hour)/controlsrc;
n     = 100;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.src = controlsrc;

cmin = (model.(itf).theta0)*(model.(itf).cmax);
vols = model.(sd).operators.vols;
% In following function, we assume that we have only one particle
computeCaverage = @(c) (sum(vols.*c)/sum(vols));
control.stopFunction = @(model, state, state0_inner) (computeCaverage(state.(sd).c) <= cmin);

schedule = struct('control', control, 'step', step); 

%% Run simulation

model.verbose = true;
[wellSols, states, report] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps', true); 

%% plotting

ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

time = cellfun(@(state) state.time, states);
cSurface = cellfun(@(state) state.(sd).cSurface, states);
phi = cellfun(@(state) state.phi, states);

figure
plot(time/hour, cSurface/(1/litre));
xlabel('time [hour]');
ylabel('Surface concentration [mol/L]');

figure
plot(time/hour, phi);
xlabel('time [hour]');
ylabel('Potential [mol/L]');

cmin = cellfun(@(state) min(state.(sd).c), states);
cmax = cellfun(@(state) max(state.(sd).c), states);

figure
hold on
plot(time/hour, cmin, 'displayname', 'cmin');
plot(time/hour, cmax, 'displayname', 'cmax');
legend show



%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
