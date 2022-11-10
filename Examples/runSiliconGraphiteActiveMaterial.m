%% run stand-alone active material model

% clear the workspace and close open figures
clear all
close all

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa

%% Setup the properties of Li-ion battery materials and cell design
jsonstruct = parseBattmoJson(fullfile('Examples', 'jsoninputs', 'silicongraphite.json'));

paramobj = SiliconGraphiteActiveMaterialInputParams(jsonstruct);

xlength = 57e-6; 
G = cartGrid(1, xlength);
G = computeGeometry(G);

gr = 'Graphite';
si = 'Silicon';

paramobj.G = G;
paramobj.(si).G = G;
paramobj.(gr).G = G;

model = SiliconGraphiteActiveMaterial(paramobj);
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

mat = gr;

sd  = 'SolidDiffusion';
itf = 'Interface';

cElectrodeInit   = (model.(mat).(itf).theta100)*(model.(mat).(itf).cmax);

mats = {gr, si};

for imat = 1 : numel(mats)
    mat = mats{imat};
    % set primary variables
    N = model.(mat).(sd).N;
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

controlsrc = 1;

total = (30*hour)/controlsrc;
n     = 100;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.src = controlsrc;

for imat = 1 : numel(mats)
    mat = mats{imat};
    cmin = (model.(mat).(itf).theta0)*(model.(mat).(itf).cmax);
    vols = model.(mat).(sd).operators.vols;
    % In following function, we assume that we have only one particle
    computeCaverage{imat} = @(c) (sum(vols.*c)/sum(vols));
end

% control.stopFunction = @(model, state, state0_inner) (feval(computeCaverage{1}, state.(mats{1}).(itf).cElectrodeSurface ) <= cmin{1} & ...
                                                      % feval(computeCaverage{2}, state.(mats{2}).(itf).cElectrodeSurface ) <= cmin{2});

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
plot(time/hour, cSurface);

figure
plot(time/hour, phi);


cmin = cellfun(@(state) min(state.(sd).c), states);
cmax = cellfun(@(state) max(state.(sd).c), states);

figure
hold on
plot(time/hour, cmin, 'displayname', 'cmin');
plot(time/hour, cmax, 'displayname', 'cmax');
legend show
