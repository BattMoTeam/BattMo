mrstModule add ad-core

an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';

filename = '/home/xavier/Matlab/Projects/battmo/ProtonicMembrane/protonicMembrane.json';
jsonstruct = fileread(filename);
jsonstruct = jsondecode(jsonstruct);

inputparams = ProtonicMembraneCellInputParams(jsonstruct);

inputparams = setupProtonicMembraneCellGrid(inputparams, jsonstruct);

% Setup model
model = ProtonicMembraneCell(inputparams);

% Setup computational graph
model = model.setupComputationalGraph();
cgt = model.computationalGraph;

% Print root variables
cgt.printRootVariables();

% Print tail variables
cgt.printTailVariables();

% plot computational graph
figure
cgt.plotComputationalGraph('nodefontsize', 14);

