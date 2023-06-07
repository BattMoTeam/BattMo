%% run stand-alone active material model

% clear the workspace and close open figures
clear all
close all

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa

%% shortcuts

gr = 'Graphite';
si = 'Silicon';

sd  = 'SolidDiffusion';
itf = 'Interface';

%% Setup the properties of Li-ion battery materials and cell design
jsonstruct = parseBattmoJson(fullfile('Examples', 'jsoninputs', 'silicongraphite.json'));

paramobj = CompositeActiveMaterialInputParams(jsonstruct);

rhoGr = paramobj.(gr).(itf).density;
rhoSi = paramobj.(si).(itf).density;

wfGr = 0.92; % weight fraction graphite
wfSi = 0.08; % weight fraction silicon

vfGr = wfGr/rhoGr;
vfSi = wfSi/rhoSi;
totV = (vfGr + vfSi);
vfGr = vfGr/totV;
vfSi = vfSi/totV;

paramobj.(gr).activeMaterialFraction      = vfGr;
paramobj.(gr).(sd).activeMaterialFraction = vfGr;
paramobj.(si).activeMaterialFraction      = vfSi;
paramobj.(si).(sd).activeMaterialFraction = vfSi;

paramobj.(gr).(itf).theta0 = 0.01;
paramobj.(si).(itf).theta0 = 0.01;
% paramobj.(si).(sd).D0 = 1e-17;
% paramobj.(si).(itf).k0 = 1e-12;
% paramobj.(gr).(sd).D0 = 1e-16;

paramobj = paramobj.validateInputParams();

xlength = 57e-6; 
G = cartGrid(1, xlength);
G = computeGeometry(G);

paramobj.G = G;
paramobj.(si).G = G;
paramobj.(gr).G = G;

model = CompositeActiveMaterial(paramobj);
model.AutoDiffBackend= AutoDiffBackend();

model.use_thermal = false;

inspectgraph = false;
if inspectgraph
    model.isRoot = true;
    cgt = ComputationalGraphTool(model);
    [g, edgelabels] = cgt.getComputationalGraph();

    figure
    % h = plot(g, 'edgelabel', edgelabels, 'nodefontsize', 10);
    h = plot(g, 'nodefontsize', 10);
    return
end

%% Setup initial state

cElectrolyte     = 5e-1*mol/litre;
phiElectrolyte   = 0;
T                = 298;

initState.cElectrolyte = cElectrolyte;
initState.phiElectrolyte = phiElectrolyte;
sd  = 'SolidDiffusion';
itf = 'Interface';

mats = {gr, si};

for imat = 1 : numel(mats)
    mat = mats{imat};
    % set primary variables
    N = model.(mat).(sd).N;
    cElectrodeInit = (model.(mat).(itf).theta0)*(model.(mat).(itf).cmax);
    initState.(mat).(sd).c        = cElectrodeInit*ones(N, 1);
    initState.(mat).(sd).cSurface = cElectrodeInit;
end

mat = gr;
% set static variable fields
initState.T = T;
initState.(mat).(itf).cElectrolyte   = cElectrolyte;
initState.(mat).(itf).phiElectrolyte = phiElectrolyte;

initState = model.dispatchTemperature(initState);
initState.(mat) = model.(mat).dispatchTemperature(initState.(mat));
initState.(mat) = model.(mat).updateConcentrations(initState.(mat));
initState.(mat).(itf) = model.(mat).(itf).updateOCP(initState.(mat).(itf));

OCP = initState.(mat).(itf).OCP;
initState.phi = OCP + phiElectrolyte;


%% setup schedule

controlsrc = -1e1;

total = (80*hour)/abs(controlsrc);
n     = 100;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.src = controlsrc;

% control.stopFunction = @(model, state, state0_inner) model.stopfunction(state, state0_inner);

schedule = struct('control', control, 'step', step); 

%% Run simulation

model.verbose = true;

nls = NonLinearSolver;
nls.errorOnFailure = false;

[wellSols, states, report] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

%% plotting

set(0, 'defaultlinelinewidth', 3);

ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

time       = cellfun(@(state) state.time, states);
cSiSurface = cellfun(@(state) state.(si).(sd).cSurface, states);
cGrSurface = cellfun(@(state) state.(gr).(sd).cSurface, states);
phi        = cellfun(@(state) state.phi, states);

figure
hold on
plot(time/hour, cSiSurface, 'displayname', 'silicon');
plot(time/hour, cGrSurface, 'displayname', 'graphite');
xlabel('time [h]')
ylabel('concentration [mol/m^3]')
legend show
title('surface concentration for the particle');

figure
hold on
plot(time/hour, cSiSurface/model.(si).(itf).cmax, 'displayname', 'silicon');
plot(time/hour, cGrSurface/model.(gr).(itf).cmax, 'displayname', 'graphite');
xlabel('time [h]')
ylabel('theta [-]')
legend show
title('surface theta for the particle');

figure
plot(time/hour, phi);
xlabel('time [h]')
ylabel('voltage [V]')
title('phi')


%%

clear c
for imat = 1 : numel(mats)
    mat = mats{imat};
    c.(mat).min = cellfun(@(state) min(state.(mat).(sd).c), states);
    c.(mat).max = cellfun(@(state) max(state.(mat).(sd).c), states);
end

matnames = {'Graphite', 'Silicon'};

for imat = 1 : numel(mats)

    mat = mats{imat};
    figure

    p1 = [time/hour, c.(mat).min];
    p2 = [time/hour, c.(mat).max];
    p2 = flip(p2, 1);
    p = [p1; p2];
    
    fill(p(:, 1), p(:, 2), 'red')
    xlabel('time [h]')
    ylabel('concentration [mol/m^3]')
    titlestr = sprintf('min and max concentration for %s', matnames{imat});
    title(titlestr)
end

