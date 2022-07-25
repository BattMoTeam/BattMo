%% run stand-alone active material model

% clear the workspace and close open figures
clear
close all

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa

%% Setup the properties of Li-ion battery materials and cell design
jsonstruct = parseBattmoJson('ParameterData/ParameterSets/Safari2009/fullmodel.json');

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

paramobj.(an).(sd).N   = 3;
paramobj.(an).(sd).np  = 1;
paramobj.(an).(sei).N  = 5;
paramobj.(an).(sei).np = 1;

xlength = 57e-6; 
G = cartGrid(1, xlength);
G = computeGeometry(G);
paramobj.(an).G = G;


paramobj.(ct).(sd).N   = 7;
paramobj.(ct).(sd).np  = 1;
paramobj.(ct).G = G;

model = SingleParticleSEI(paramobj);

model.(ctrl).Imax = 1e-3;

dograph = false;

if dograph
    model = model.registerVarAndPropfuncNames();
    [g, edgelabels] = setupGraph(model);
    cgf = ComputationalGraphFilter(g);
    % cgf.includeNodeNames = 'Anode.SolidDiffusion.cSurface';
    % cgf.includeNodeNames = 'cInterface';
    % g = cgf.setupGraph('oneParentOnly', true);
    % gg = cgf.setupDescendantGraph();
    % g = cgf.setupGraph('oneParentOnly', true);    
    figure
    h = plot(gg); 
end

%% Setup initial state

NanodeSd  = model.(an).(sd).N;
NanodeSEI = model.(an).(sei).N;
cAnodeMax = model.(an).(itf).cmax;

NcathodeSd  = model.(ct).(sd).N;
cCathodeMax = model.(ct).(itf).cmax;

x0 = 0.75;                     % Initial stochiometry from Safari
cAnode = x0*cAnodeMax;
y0 = 0.5;                     % Initial stochiometry from Safari
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
initState.(an).(sei).delta      = 0;
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


%% setup schedule

total = 1*hour;
n     = 100;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

switch model.Control.controlPolicy
  case 'IEswitch'
    tup = 0.1; % rampup value for the current function, see rampupSwitchControl
    srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                model.Control.Imax, ...
                                                model.Control.lowerCutoffVoltage);
    % we setup the control by assigning a source and stop function.
    control = struct('src', srcfunc, 'IEswitch', true);
  case 'CCCV'
    control = struct('CCCV', true);
  otherwise
    error('control policy not recognized');
end

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step);


%% Run simulation

model.verbose = true;

dopack = false;
if dopack
    dataFolder = 'BattMo';
    problem = packSimulationProblem(initState, model, schedule, dataFolder, 'Name', 'temp');
    problem.SimulatorSetup.OutputMinisteps = true; 
    simulatePackedProblem(problem);
    [globvars, states, report] = getPackedSimulatorOutput(problem);
else
    [wellSols, states, report] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps', true); 
end

%% Plotting

% concentration evolution in particle

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

%%

figure
d = cellfun(@(state) state.(sei).delta, states);
t = cellfun(@(state) state.time, states);
plot(t, d);
xlabel('time / [s]');
ylabel('sie width / [m]');






