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
cgt = model.computationalGraph;

% Print root variables
cgt.printRootVariables();

% Print tail variables
cgt.printTailVariables();

% plot computational graph
[g, edgelabels] = cgt.getComputationalGraph();

figure
h = plot(g, 'nodefontsize', 14);
