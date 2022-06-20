mrstModule add ad-core

an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';

paramobj = ProtonicMembraneCellInputParams([]);

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

model = ProtonicMembraneCell(paramobj);

model = model.registerVarAndPropfuncNames();

[g, edgelabels] = setupGraph(model);

cgf = ComputationalGraphFilter(g);
cgf.includeNodeNames = [];

g = cgf.setupGraph();

figure
h = plot(g, 'nodefontsize', 18);