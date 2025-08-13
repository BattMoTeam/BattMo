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
cgti = model.cgti;

% Print root variables
cgti.printRootVariables();

% Print tail variables
cgti.printTailVariables();

% plot computational graph
figure
cgti.plotComputationalGraph('nodefontsize', 14);

