%% run stand-alone active material model

% clear the workspace and close open figures
clear all
close all

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa

%% Setup the properties of Li-ion battery materials and cell design
jsonstruct = parseBattmoJson(fullfile('ParameterData','ParameterSets','Safari2009','anode_sei.json'));

% Some shorthands used for the sub-models
ne    = 'NegativeElectrode';
pe    = 'PositiveElectrode';
am    = 'ActiveMaterial';
sd    = 'SolidDiffusion';
itf   = 'Interface';
sei   = 'SolidElectrodeInterface';
sr    = 'SideReaction';
elyte = 'Electrolyte';

paramobj = SEIActiveMaterialInputParams(jsonstruct);

paramobj.(sd).N  = 10;
paramobj.(sei).N = 10;

paramobj.standAlone = true;

paramobj = paramobj.validateInputParams();

model = SEIActiveMaterial(paramobj);
% model = ActiveMaterial(paramobj);
% model = SolidElectrodeInterface(paramobj.(sei));
% model = SideReaction(paramobj.(sr));
% model = Interface(paramobj.(itf));
% model = FullSolidDiffusionModel(paramobj.(sd));

doplotgraph = false;
if doplotgraph
    cgt = model.computationalGraph;
    if isempty(cgt)
        model = model.setupComputationalGraph();
        cgt = model.computationalGraph;
    end
    cgt.plotComputationalGraph()
    return
end

%% Setup initial state

Nsd  = model.(sd).N;
Nsei = model.(sei).N;

cElectrodeInit   = 0.75*model.(itf).saturationConcentration;
phiElectrodeInit = 0;
cElectrolyte     = 5e-1*mol/litre;
T                = 298.15; % reference temperature in Interface.

epsiSEI     = 0.05;            % From Safari
cECsolution = 4.541*mol/litre; % From Safari
cECexternal = epsiSEI*cECsolution;

% compute OCP and  phiElectrolyte
initState.T = T;
initState.(sd).cSurface = cElectrodeInit;
initState = model.evalVarName(initState, {itf, 'OCP'});

OCP = initState.(itf).OCP;
phiElectrolyte = phiElectrodeInit - OCP;

% set primary variables
initState.E                = phiElectrodeInit;
initState.I                = 0;
initState.(sd).c           = cElectrodeInit*ones(Nsd, 1);
initState.(sei).c          = cECexternal*ones(Nsei, 1);
initState.(sei).cInterface = cECexternal;
initState.(sei).delta      = 5*nano*meter;
initState.R                = 0;

% set static variable fields
initState.T = T;
initState.(itf).cElectrolyte   = cElectrolyte;
initState.(itf).phiElectrolyte = phiElectrolyte;
initState.(sr).phiElectrolyte  = phiElectrolyte;
initState.(sei).cExternal      = cECexternal;


%% setup schedule

% Reference rate which roughly corresponds to 1 hour for the data of this example
Iref = 1.3e-4*ampere/(1*centi*meter)^2;

Imax = 1e1*Iref;

total = 1*hour*(Iref/Imax);
n     = 100;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

tup = dt; % rampup value for the current function, see rampupSwitchControl
srcfunc = @(time) rampupControl(time, tup, Imax);

cmin = (model.(itf).guestStoichiometry0)*(model.(itf).saturationConcentration);
control.stopFunction = @(model, state, state0_inner) (state.(sd).cSurface <= cmin);
control.src = srcfunc;

schedule = struct('control', control, 'step', step); 

%% setup non-linear solver

nls = NonLinearSolver();
nls.errorOnFailure = false;

model.nonlinearTolerance = 1e-5;

%% Run simulation

model.verbose = true;
[wellSols, states, report] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

%% Plotting

set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);
set(0, 'defaultlinelinewidth', 3);

ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

time = cellfun(@(state) state.time, states);

cSurface = cellfun(@(state) state.(sd).cSurface, states);
figure
plot(time/hour, cSurface/(1/litre));
xlabel('time / h');
ylabel('Surface concentration / mol/L');
title('Surface concentration');

E = cellfun(@(state) state.E, states);
figure
plot(time/hour, E);
xlabel('time / h');
ylabel('Potential / V');
title('Potential');

I = cellfun(@(state) state.I, states);
figure
plot(time/hour, I);
xlabel('time / h');
ylabel('Current density / A/m^2');
title('Current density');


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

delta = cellfun(@(state) state.(sei).delta, states);
figure
plot(time/hour, delta/(nano*meter));
xlabel('time [hour]');
ylabel('thickness / nm');
title('SEI thickness')

c = states{end}.(sd).c;
r = linspace(0, model.(sd).particleRadius, model.(sd).N);

figure
plot(r, c/(mol/litre));
xlabel('radius / m')
ylabel('concentration / mol/L')
title('Particle concentration profile (last time step)')

r = states{end}.(sei).delta;
r = linspace(0, r, model.(sei).N);
c = states{end}.(sei).c;

figure
plot(r/(nano*meter), c/(mol/litre));
xlabel('x / mm')
ylabel('concentration / mol/L');
title('Concentration profile in SEI layer (last time step)');






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
