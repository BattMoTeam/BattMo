mrstModule add ad-core

an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';

filename = '/home/xavier/Matlab/Projects/battmo/ProtonicMembrane/protonicMembrane.json';
jsonstruct = fileread(filename);
jsonstruct = jsondecode(jsonstruct);

paramobj = ProtonicMembraneCellInputParams(jsonstruct);

G = cartGrid(1, 10);
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

% Setup computational graph
model = model.setupComputationalGraph();
cgt = model.computationalGraph;

% Print root variables
cgt.printRootVariables();

% Print tail variables
cgt.printTailVariables();

% plot computational graph
[g, edgelabels] = cgt.getComputationalGraph();

figure
h = plot(g, 'nodefontsize', 14);
