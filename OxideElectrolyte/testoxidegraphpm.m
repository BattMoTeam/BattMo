mrstModule add ad-core

an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';

filename = '/home/xavier/Matlab/Projects/battmo/OxideMembrane/protonicMembrane.json';
jsonstruct = fileread(filename);
jsonstruct = jsondecode(jsonstruct);

paramobj = OxideMembraneCellInputParams(jsonstruct);

paramobj = setupOxideMembraneCellGrid(paramobj, jsonstruct);

% Setup model
model = OxideMembraneCell(paramobj);

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
