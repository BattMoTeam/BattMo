mrstModule add ad-core

an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';

filename = '/home/xavier/Matlab/Projects/battmo/OxideElectrolyte/oxidemembrane.json';
jsonstruct = fileread(filename);
jsonstruct = jsondecode(jsonstruct);

inputparams = OxideMembraneCellInputParams(jsonstruct);

inputparams = setupProtonicMembraneCellGrid(inputparams, jsonstruct);

% Setup model
model = OxideMembraneCell(inputparams);

% Setup computational graph
model = model.setupComputationalGraph();
cgit = model.computationalGraph;

% Print root variables
cgit.printRootVariables();

% Print tail variables
cgit.printTailVariables();

% plot computational graph
[g, edgelabels] = cgit.getComputationalGraph();

figure
h = plot(g, 'nodefontsize', 14);
