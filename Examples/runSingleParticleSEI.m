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

paramobj = SingleParticleSEIInputParams(jsonstruct);

paramobj.(an).(sd).N   = 10;
paramobj.(an).(sd).np  = 1;
paramobj.(an).(sei).N  = 10;
paramobj.(an).(sei).np = 1;

xlength = 57e-6; 
G = cartGrid(1, xlength);
G = computeGeometry(G);
paramobj.(an).G = G;


paramobj.(ct).(sd).N   = 10;
paramobj.(ct).(sd).np  = 1;
paramobj.(ct).G = G;

model = SingleParticleSEI(paramobj);

dograph = true;

if dograph
    model = model.registerVarAndPropfuncNames();
    [g, edgelabels] = setupGraph(model);
    % cgf = ComputationalGraphFilter(g);
    % cgf.includeNodeNames = 'external';
    % g = cgf.setupGraph();
    
    figure
    h = plot(g); 
end

return

%% Setup initial state

Nsd = model.(sd).N;
Nsei = model.(sei).N;

cElectrodeInit   = 40*mol/litre;
phiElectrodeInit = 3.5;
cElectrolyte     = 5e-1*mol/litre;
T                = 298.15; % reference temperature in Interface.
cExternal        = 4.541*mol/litre;
% compute OCP and  phiElectrolyte
clear state
state.cElectrodeSurface = cElectrodeInit;
state.T = T;
state = model.(itf).updateOCP(state);
OCP = state.OCP;
phiElectrolyte = phiElectrodeInit - OCP;

% set primary variables
initState.phi              = phiElectrodeInit;
initState.(sd).c           = cElectrodeInit*ones(Nsd, 1);
initState.(sd).cSurface    = cElectrodeInit;
initState.(sei).c          = cElectrolyte*ones(Nsei, 1);
initState.(sei).cInterface = cElectrolyte;
initState.(sei).delta      = 0;
initState.R                = 0;

% set static variable fields
initState.T = T;
initState.(itf).cElectrolyte   = cElectrolyte;
initState.(itf).phiElectrolyte = phiElectrolyte;
initState.(sr).phiElectrolyte  = phiElectrolyte;
initState.(sei).cExternal      = cExternal;


%% setup schedule

total = 1*hour;
n     = 100;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.src = 1;

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






