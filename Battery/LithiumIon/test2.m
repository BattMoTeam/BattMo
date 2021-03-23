% 1D test

clear all
close all

%% add MRST module
mrstModule add ad-core multimodel

G = cartGrid(10, 10);
G = computeGeometry(G);

nc = G.cells.num;

cells = (1 : G.cells.num)';

model = Battery();

elyte = model.getAssocModel('elyte');
ne    = model.getAssocModel('ne');
pe    = model.getAssocModel('pe');
ccne  = model.getAssocModel('ccne');
ccpe  = model.getAssocModel('ccpe');

n = 0;

nc = elyte.G.cells.num;
state.elyte_phi = ones(nc, 1);
state.elyte_cs = cell(1, 2);
state.elyte_cs{1} = ones(nc, 1);
n = n + nc;

nc = ne.G.cells.num;
state.ne_am_phi = ones(nc, 1);
state.ne_am_Li = ones(nc, 1);
n = n + nc;

nc = pe.G.cells.num;
state.pe_am_phi = ones(nc, 1);
state.pe_am_Li = ones(nc, 1);
n = n + nc;

nc = ccne.G.cells.num;
state.ccne_phi = ones(nc, 1);
n = n + nc;

nc = ccpe.G.cells.num;
state.ccpe_phi = ones(nc, 1);
state.ccpe_E = 1;
n = n + nc;

model = model.setupFV(state);
fv = model.fv;

y = NaN(nc, 1);

pvarnames = model.getModelPrimaryVarNames();

for ind = 1 : numel(pvarnames)
    varname = pvarnames{ind};
    val = model.getProp(state, varname);
    y(fv.getSlot(varname)) = val;
end

    


