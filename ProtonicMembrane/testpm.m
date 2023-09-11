clear all

mrstModule add ad-core

an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';

filename = '/home/xavier/Matlab/Projects/battmo/ProtonicMembrane/protonicMembrane.json';
jsonstruct = fileread(filename);
jsonstruct = jsondecode(jsonstruct);

paramobj = ProtonicMembraneCellInputParams(jsonstruct);

G = cartGrid(10, 1);
G = computeGeometry(G);

paramobj.(elyte).G = G;

couplingTerms = {};

coupterm = couplingTerm('Anode-Electrolyte', {an, elyte});
coupterm.couplingcells = [1, 1];
coupterm.couplingfaces = [1, 1];
couplingTerms{end + 1} = coupterm;

coupterm = couplingTerm('Cathode-Electrolyte', {an, elyte});
coupterm.couplingcells = [1, G.cells.num];
coupterm.couplingfaces = [1, G.faces.num];
couplingTerms{end + 1} = coupterm;

paramobj.couplingTerms = couplingTerms;

% Setup model
model = ProtonicMembraneCell(paramobj);

model.verbose = true;

% Setup initial state
state0 = model.setupInitialState();

% Setup fake schedule (not needed here)
step = struct('val', 1, 'control', 1);
control.dummy = [];
schedule = struct('control', control, 'step', step); 

[~, states, report] = simulateScheduleAD(state0, model, schedule); 



