mrstModule add ad-core

an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';

filename = '/home/xavier/Matlab/Projects/battmo/ProtonicMembrane/protonicMembrane.json';
jsonstruct = fileread(filename);
jsonstruct = jsondecode(jsonstruct);

inputparams = ProtonicMembraneInputParams(jsonstruct);

inputparams = setupProtonicMembraneGrid(inputparams, jsonstruct);

% Setup model
model = ProtonicMembrane(inputparams);

% Setup computational graph
cgit = model.cgit;

% Print root variables
cgit.printRootVariables();

% Print tail variables
cgit.printTailVariables();

% plot computational graph
figure
cgit.plotComputationalGraph('nodefontsize', 14);

