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

%% script to order the property function 
% require bgl

% addpath(genpath('/home/xavier/Programs/matlab_bgl/'));

[g, edgelabels] = setupGraph(model);

nn = g.numnodes;
nodenames = g.Nodes.Variables;


[c, ia, ic] = unique(edgelabels.ts);

propnames = c;
propindex = edgelabels.ps(ia);
propindex = cell2mat(propindex);

A = adjacency(g);

p = topological_order(A);

nodenames = nodenames(p);

[lia, locb] = ismember(nodenames, propnames);
propindex = propindex(locb(lia));

for ind = 1 : numel(propindex)
    iprop = propindex(ind);
    propfunction = model.propertyFunctionList{iprop};
    fn = propfunction.fn;
    fnname = func2str(fn);
    fprintf('%s\n', fnname);
end
