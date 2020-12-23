% 1D test

clear all
close all

%% add MRST module
mrstModule add ad-core

G = cartGrid(10, 10);
G = computeGeometry(G);

nc = G.cells.num;

cells = (1 : G.cells.num)';

state.T = ones(nc, 1);
state.SOC = ones(nc, 1);
state.elyte_c_Li = ones(nc, 1);
state.c_Li = ones(nc, 1);

model = orgLiPF6('elyte', G, cells);
% model = BatteryModel();

state = model.initializeState(state);
            
