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
    cgt.plotComputationalGraph()
    return
end

% cgt.printTailVariables
% return

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

Imax = Iref;

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
[wellSols, states, report] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

%% Plotting

ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

time     = cellfun(@(state) state.time, states);

cSurface = cellfun(@(state) state.(sd).cSurface, states);

figure
plot(time/hour, cSurface/(1/litre));
xlabel('time [hour]');
ylabel('Surface concentration [mol/L]');

E        = cellfun(@(state) state.E, states);

figure
plot(time/hour, E);
xlabel('time [hour]');
ylabel('Potential / Voltage');

return

%%


ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

doplotconcs = false;
% concentration evolution in particle
if doplotconcs
    figure
    xr = (model.(sd).rp/model.(sd).N) * (1 : model.(sd).N)';
    for ind = 1 : numel(states)
        state = states{ind};
        cla
        plot(xr, state.(sd).c)
        xlabel('r / [m]')
        ylabel('concentration')
        title(sprintf('time : %g s', state.time));
        pause(0.1)
    end
end

%%

figure
d = cellfun(@(state) state.(sei).delta, states);
t = cellfun(@(state) state.time, states);
plot(t, d);
xlabel('time / [s]');
ylabel('sie width / [m]');


figure
E = cellfun(@(state) state.E, states);
plot(t, E);
xlabel('time / [s]');
ylabel('voltage / [V]');








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
