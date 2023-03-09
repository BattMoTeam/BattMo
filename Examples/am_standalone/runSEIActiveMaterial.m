%% run stand-alone active material model

% clear the workspace and close open figures
clear
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

%% We setup the battery geometry ("bare" battery with no current collector). We use that to recover the parameters for the active material of the positive electrode, which is instantiate later
paramobj = SEIActiveMaterialInputParams(jsonstruct);

paramobj.externalCouplingTerm = [];
paramobj.(sd).N   = 10;
paramobj.(sd).np  = 1;
paramobj.(sei).N  = 10;
paramobj.(sei).np = 1;
xlength = 57e-6; 
G = cartGrid(1, xlength);
G = computeGeometry(G);
paramobj.G = G;

model = SEIActiveMaterial(paramobj);
model.standAlone = true;

cgt = ComputationalGraphTool(model);
% cgt.printTailVariables
% return

%% Setup initial state

Nsd = model.(sd).N;
Nsei = model.(sei).N;

cElectrodeInit   = 0.75*model.(itf).cmax;
phiElectrodeInit = 0;
cElectrolyte     = 5e-1*mol/litre;
T                = 298.15; % reference temperature in Interface.

epsiSEI     = 0.05;                % From Safari
cECsolution = 4.541*mol/litre; % From Safari
cECexternal = epsiSEI*cECsolution;

% compute OCP and  phiElectrolyte
clear state
state.cElectrodeSurface = cElectrodeInit;
state.T = T;
state = model.(itf).updateOCP(state);
OCP = state.OCP;
phiElectrolyte = phiElectrodeInit - OCP;

% set primary variables
initState.E                = phiElectrodeInit;
initState.(sd).c           = cElectrodeInit*ones(Nsd, 1);
initState.(sd).cSurface    = cElectrodeInit;
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

total = 1*hour;
n     = 100;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.src = 2e-5;

schedule = struct('control', control, 'step', step); 

%% Run simulation

model.verbose = true;

nls = NonLinearSolver;
nls.errorOnFailure = false;

dopack = false;
if dopack
    dataFolder = 'BattMo';
    problem = packSimulationProblem(initState, model, schedule, dataFolder, 'Name', 'temp');
    problem.SimulatorSetup.OutputMinisteps = true; 
    simulatePackedProblem(problem);
    [globvars, states, report] = getPackedSimulatorOutput(problem);
else
    [wellSols, states, report] = simulateScheduleAD(initState, model, schedule, ...
                                                    'OutputMinisteps', true, ...
                                                    'NonLinearSolver', nls); 
end

%% Plotting


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





