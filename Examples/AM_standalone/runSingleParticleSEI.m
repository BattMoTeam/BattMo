%% run stand-alone active material model

% clear the workspace and close open figures
clear
close all

%% Import the required modules from MRST
% load MRST modules

mrstModule add ad-core mrst-gui mpfa

%% Setup the properties of Li-ion battery materials and cell design
jsonstruct = parseBattmoJson(fullfile('ParameterData', 'ParameterSets', 'Safari2009', 'fullmodel.json'));

% Some shorthands used for the sub-models
an    = 'Anode';
ct    = 'Cathode';
sd    = 'SolidDiffusion';
itf   = 'Interface';
sei   = 'SolidElectrodeInterface';
sr    = 'SideReaction';
elyte = 'Electrolyte';
ctrl  = 'Control';

inputparams = SingleParticleSEIInputParams(jsonstruct);

controlPolicy = 'CCCV';
switch controlPolicy
  case 'CCCV'
    inputparams.Control.controlPolicy = 'CCCV';
  case 'CV'
    inputparams_control = CvControlModelInputParams([]);
    inputparams_control.inputVoltage = 4.5;
    inputparams.Control = inputparams_control;
  otherwise
    error('control policy not recognized')
end

inputparams.(an).(sd).N   = 10;
inputparams.(an).(sd).np  = 1;
inputparams.(an).(sei).N  = 10;
inputparams.(an).(sei).np = 1;

xlength = 57e-6;
G = cartGrid(1, xlength);
inputparams.(an).G = Grid(G);

inputparams.(ct).(sd).N  = 10;
inputparams.(ct).(sd).np = 1;
inputparams.(ct).G = Grid(G);

model = SingleParticleSEI(inputparams);
model = model.equipModelForComputation();

% Normally in CCCV control, the current value is computed from DRate but, here, we set it up directly
if strcmp(model.(ctrl).controlPolicy, 'CCCV')
    model.(ctrl).ImaxCharge    = 1.1;
    model.(ctrl).ImaxDischarge = 1.8;
end

dograph = false;

if dograph
    cgti = ComputationalGraphInteractiveTool(model);
    % cgti.includeNodeNames = 'Anode.SolidDiffusion.cSurface';
    cgti.includeNodeNames = 'Control.I$';
    % cgti.includeNodeNames = 'SideReaction.R';
    % [g, edgelabels] = cgti.setupGraph('oneParentOnly', true, 'type', 'ascendant');
    % [g, edgelabels] = cgti.setupGraph('oneParentOnly', true, 'type', 'descendant');
    % [g, edgelabels] = cgti.setupGraph();
    [g, edgelabels] = cgti.setupDescendantGraph();
    figure
    h = plot(g, 'edgelabel', edgelabels,'edgefontsize', 15, 'nodefontsize', 12);
    return
end


%% Setup initial state

NanodeSd  = model.(an).(sd).N;
NanodeSEI = model.(an).(sei).N;
cAnodeMax = model.(an).(itf).saturationConcentration;

NcathodeSd  = model.(ct).(sd).N;
cCathodeMax = model.(ct).(itf).saturationConcentration;

x0 = 0.75; % Initial stochiometry from Safari
cAnode = x0*cAnodeMax;
y0 = 0.5; % Initial stochiometry from Safari
cCathode = y0*cCathodeMax;

phiAnode = 0;
T = 298.15;                    % NOTE : should match reference temperature in Interface in order to reproduce Safari result.
epsiSEI = 0.05;                % From Safari
cECsolution = 4.541*mol/litre; % From Safari
cECexternal = epsiSEI*cECsolution;

% compute OCP at anode and  phiElectrolyte
clear an_itf_state
an_itf_state.cElectrodeSurface = cAnode;
an_itf_state.T = T;
an_itf_state = model.(an).(itf).updateOCP(an_itf_state);
OCP = an_itf_state.OCP;
phiElectrolyte = phiAnode - OCP;

% compute OCP at cathode and  phiCathode
clear ct_itf_state
ct_itf_state.cElectrodeSurface = cCathode;
ct_itf_state.T = T;
ct_itf_state = model.(ct).(itf).updateOCP(ct_itf_state);
OCP = ct_itf_state.OCP;
phiCathode = phiElectrolyte + OCP;

% set primary variables
initState.(an).(sd).c           = cAnode*ones(NanodeSd, 1);
initState.(an).(sd).cSurface    = cAnode;
initState.(an).(sei).c          = cECexternal*ones(NanodeSEI, 1);
initState.(an).(sei).cInterface = cECexternal;
initState.(an).(sei).delta      = 5*nano*meter;
initState.(an).R                = 0;

initState.(ct).(sd).c           = cCathode*ones(NcathodeSd, 1);
initState.(ct).(sd).cSurface    = cCathode;

initState.(elyte).phi = phiElectrolyte;

initState.(ctrl).E = phiCathode;
initState.(ctrl).I = 0;

% set static variable fields
initState.(an).phi                         = 0;
initState.T                                = T;
initState.(an).(sei).cExternal             = cECexternal;
initState.(ct).(itf).externalPotentialDrop = 0;

switch model.(ctrl).controlPolicy
  case 'CCCV'
    initState.(ctrl).E = 0;
    initState.(ctrl).ctrlType     = 'CV_charge2';
  case 'CV'
    % ok. nothing to do
  otherwise
    error('control policy not recognized');
end

initstate.time = 0;

%% setup schedule

total = 800*hour;
n     = 50*800;
dt    = total/n;

n  = 400;
dt = dt*ones(n, 1);


if strcmp(model.Control.controlPolicy, 'CV')
    dt = rampupTimesteps(800*hour, 1*hour, 4);
end

step  = struct('val', dt, 'control', ones(size(dt, 1), 1));

switch model.Control.controlPolicy
  case 'CCCV'
    control = struct('CCCV', true);
  case 'CV'
    control = struct('CV', true);
  otherwise
    error('control policy not recognized');
end

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step);


%% Run simulation

model.verbose = true;

% Setup nonlinear solver
nls = NonLinearSolver();
nls.maxTimestepCuts = 10;
% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 100;
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false;

dopack = false;
if dopack
    dataFolder = 'BattMo';
    problem = packSimulationProblem(initState, model, schedule, dataFolder, 'Name', 'safari3', 'NonLinearSolver', nls);
    problem.SimulatorSetup.OutputMinisteps = true;
    simulatePackedProblem(problem);
    [globvars, states, report] = getPackedSimulatorOutput(problem);
else
    [~, states, report] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);
end

%% Process output and recover the output voltage and current from the output states.

ind = cellfun(@(x) not(isempty(x)), states);
states = states(ind);
E = cellfun(@(x) x.Control.E, states);
I = cellfun(@(x) x.Control.I, states);
% [SOCN, SOCP] =  cellfun(@(x) model.calculateSOC(x), states);
time = cellfun(@(x) x.time, states);

figure
plot(time/hour, E);
xlabel('time / [h]');
ylabel('Voltage / [V]');

figure
d = cellfun(@(state) state.(an).(sei).delta, states);
t = cellfun(@(state) state.time, states);
plot(t/hour, d/(micro*meter));
xlabel('time / [h]');
ylabel('sie width / [\mu m]');



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
