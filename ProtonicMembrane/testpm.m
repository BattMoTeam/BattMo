clear all
close all

mrstModule add ad-core

an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';
ctrl  = 'Control';

filename = '/home/xavier/Matlab/Projects/battmo/ProtonicMembrane/protonicMembrane.json';
jsonstruct = fileread(filename);
jsonstruct = jsondecode(jsonstruct);

paramobj = ProtonicMembraneCellInputParams(jsonstruct);

paramobj = setupProtonicMembraneCellGrid(paramobj, jsonstruct);

% Setup model
model = ProtonicMembraneCell(paramobj);

model = model.setupComputationalGraph();

% compute and get computationalGraph (just used for postprocessing)
model = model.validateModel();
cgt   = model.computationalGraph;

model.verbose = true;

% Setup initial state
state0 = model.setupInitialState();

% Setup schedule

T = 1; % This is not a real time scale, as all the model deals with equilibrium
N = 10;
dt = T/N;

step.val = dt*ones(N, 1);
step.control = ones(N, 1);

Imax = 1e-1;

control.src = @(time) time/T*Imax;

schedule = struct('control', control, 'step', step); 

nls = NonLinearSolver();
nls.maxIterations = 100;
nls.errorOnFailure = false;

[~, states, report] = simulateScheduleAD(state0, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

state = states{1};
state = model.addVariables(state, control.src);


