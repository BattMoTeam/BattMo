clear all
close all

mrstModule add ad-core

an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';

filename = '/home/xavier/Matlab/Projects/battmo/ProtonicMembrane/protonicMembrane.json';
jsonstruct = fileread(filename);
jsonstruct = jsondecode(jsonstruct);

paramobj = ProtonicMembraneCellInputParams(jsonstruct);

paramobj = setupProtonicMembraneCellGrid(paramobj, jsonstruct);

% Setup model
model = ProtonicMembraneCell(paramobj);

model = model.setupComputationalGraph();
model = model.validateModel();

model.verbose = true;

% Setup initial state
state0 = model.setupInitialState();

% Setup schedule

T = 1; % This is not a real time scale, as all the model deals with equilibrium
N = 100;
dt = T/N;

step.val = dt*ones(N, 1);
step.control = ones(N, 1);

Imax = 1e-2;

control.src = @(time) time/T*Imax;

schedule = struct('control', control, 'step', step); 

nls = NonLinearSolver();
nls.maxIterations = 100;

[~, states, report] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls); 

state = states{end};
state = model.addVariables(state);


