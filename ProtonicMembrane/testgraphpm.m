mrstModule add ad-core

paramobj = ProtonicMembraneCellInputParams([]);

model = ProtonicMembraneCell(paramobj);

model = model.registerVarAndPropfuncNames();

[g, edgelabels] = setupGraph(model);

cgf = ComputationalGraphFilter(g);
cgf.includeNodeNames = [];

g = cgf.setupGraph();

figure
h = plot(g, 'nodefontsize', 18);