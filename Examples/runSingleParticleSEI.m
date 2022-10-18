%% run stand-alone active material model

% clear the workspace and close open figures
clear all
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

paramobj = SingleParticleSEIInputParams(jsonstruct);

controlPolicy = 'CCCV';
switch controlPolicy 
  case 'CCCV'
    paramobj.Control.controlPolicy = 'CCCV';
  case 'CV'
    paramobj_control = CvControlModelInputParams([]);
    paramobj_control.inputVoltage = 4.5;
    paramobj_control.CRate = 0.5; % not used in this case
    paramobj.Control = paramobj_control;
  otherwise
    error('control policy not recognized')
end

paramobj.(an).(sd).N   = 10;
paramobj.(an).(sd).np  = 1;
paramobj.(an).(sei).N  = 10;
paramobj.(an).(sei).np = 1;

xlength = 57e-6; 
G = cartGrid(1, xlength);
G = computeGeometry(G);
paramobj.(an).G = G;

paramobj.(ct).(sd).N  = 10;
paramobj.(ct).(sd).np = 1;
paramobj.(ct).G = G;

model = SingleParticleSEI(paramobj);

% Normally in CCCV control, the current value is computed from CRate but, here, we set it up directly
if strcmp(model.(ctrl).controlPolicy, 'CCCV')
    model.(ctrl).Imax = 1.8;
end

dograph = false;

if dograph
    cgt = ComputationalGraphTool(model);
    % cgt.includeNodeNames = 'Anode.SolidDiffusion.cSurface';
    cgt.includeNodeNames = 'Control.I$';
    % cgt.includeNodeNames = 'SideReaction.R';
    % [g, edgelabels] = cgt.setupGraph('oneParentOnly', true, 'type', 'ascendant');
    % [g, edgelabels] = cgt.setupGraph('oneParentOnly', true, 'type', 'descendant');
    % [g, edgelabels] = cgt.setupGraph();
    [g, edgelabels] = cgt.setupDescendantGraph();
    figure
    h = plot(g, 'edgelabel', edgelabels,'edgefontsize', 15, 'nodefontsize', 12);
    return
end


%% Setup initial state

NanodeSd  = model.(an).(sd).N;
NanodeSEI = model.(an).(sei).N;
cAnodeMax = model.(an).(itf).cmax;

NcathodeSd  = model.(ct).(sd).N;
cCathodeMax = model.(ct).(itf).cmax;

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
    initState.(ctrl).ctrlType = 'CV_charge2';
    initState.(ctrl).nextCtrlType = 'CV_charge2';
  case 'CV'
    % ok. nothing to do
  otherwise
    error('control policy not recognized');
end

%% setup schedule

total = 800*hour;
n     = 50*800;
dt    = total/n;

n = 400;
dt    = dt*ones(n, 1);


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
    problem = packSimulationProblem(initState, model, schedule, dataFolder, 'Name', 'safari3', 'NonlinearSolver', nls);
    problem.SimulatorSetup.OutputMinisteps = true; 
    simulatePackedProblem(problem);
    [globvars, states, report] = getPackedSimulatorOutput(problem);
else
    [wellSols, states, report] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps', true, 'NonlinearSolver', nls); 
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
