clear all

mrstModule add ad-core mrst-gui

gen = GasSupplyPEMgridGenerator2D();

gen.nxCell     = 10;
gen.nxGasSupply = 10;

gen.lxCell     = 22*micro*meter;
gen.lxGasSupply = 1*centi*meter;

gen.ny = 10;
gen.ly = 1;

paramobj = ProtonicMembraneCellWithGasSupplyInputParams([]);

paramobj = gen.updateInputParams(paramobj);


return
%%  setup gas supply

filename = '/home/xavier/Matlab/Projects/battmo/ProtonicMembrane/gaslayersupply.json';
jsonstruct = fileread(filename);
jsonstruct = jsondecode(jsonstruct);

paramobj = ProtonicMembraneGasSupplyInputParams(jsonstruct);

gen = GasSupplyGridGenerator2D();

gen.nx = 100;
gen.ny = 70;
gen.lx = 10;
gen.ly = 10;

paramobj = gen.updateInputParams(paramobj);

paramobjGasSupply = paramobj;

%% setup membrane

an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';
ctrl  = 'Control';

filename = '/home/xavier/Matlab/Projects/battmo/ProtonicMembrane/protonicMembrane.json';
jsonstruct = parseBattmoJson(filename);

jsonstruct.(elyte).N = 10000;

paramobj = ProtonicMembraneCellInputParams(jsonstruct);

paramobj = setupProtonicMembraneCellGrid(paramobj, jsonstruct);

paramobjCell = paramobj;

%%

clear paramobj

paramobj = ProtonicMembraneCellWithGasSupplyInputParams([]);

paramobj.Cell      = paramobjCell;
paramobj.GasSupply = paramobjGasSupply;

model = ProtonicMembraneCellWithGasSupply(paramobj);

model = model.setupComputationalGraph();
cgt = model.computationalGraph;



