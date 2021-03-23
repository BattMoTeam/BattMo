% 1D test

clear all
close all

%% add MRST module
mrstModule add ad-core multimodel

G = cartGrid(10, 10);
G = computeGeometry(G);

nc = G.cells.num;

cells = (1 : G.cells.num)';

modelcase = 'battery';

switch modelcase
  case 'elyte'
    state.T = ones(nc, 1);
    state.SOC = ones(nc, 1);
    state.phi = ones(nc, 1);
    cs = cell(1, 2);
    cs{1} = ones(nc, 1);
    cs{2} = [];
    state.cs = cs;
    model = orgLiPF6('elyte', G, cells);
  case 'graphite'
    state.T = ones(nc, 1);
    state.SOC = ones(nc, 1);
    state.phielyte = ones(nc, 1);
    model = graphiteElectrode('ne', G, cells);
  case 'nmc111'
    state.T = ones(nc, 1);
    state.SOC = ones(nc, 1);
    state.phielyte = ones(nc, 1);
    model = nmc111Electrode('pe', G, cells);
  case 'currentcollector'
    state.T = ones(nc, 1);
    state.SOC = ones(nc, 1);
    model = currentCollector('ccpe', G, cells);
  case 'battery'
    
    model = Battery();
    adminmodel = AdminModel();
    model = model.setupAdminModel(adminmodel);

    nc = model.G.cells.num;
    state.T = ones(nc, 1);
    state.SOC = ones(nc, 1);
    
    % initialize elyte
    elyte = model.getAssocModel('elyte');
    nc = elyte.G.cells.num;
    state.elyte_phi = ones(nc, 1);
    cs = cell(1, 2);
    cs{1} = ones(nc, 1);
    cs{2} = [];
    state.elyte_cs = cs;
    
end

state = model.initializeState(state);

return

switch modelcase
  case {'elyte', 'currentCollector'}
    state = model.updateProp(state, 'j');
  case {'graphite', 'nmc111'}
    state = model.updateProp(state, 'R');
    state = model.updateProp(state, 'j');
end

