
clear all
close all

%% Add MRST module
mrstModule add ad-core multimodel

modelcase = 'battery';

switch modelcase
  case 'elyte'
    G = cartGrid([10, 10]);
    G = computeGeometry(G);
    model = orgLiPF6('elyte', G, (1 : G.cells.num));
  case 'battery'
    model = BatteryModel();
    model.J = 0.1;
  otherwise
    error('modelcase not recognized');
end

[g, edgelabels] = setupGraph(model);

figure
h = plot(g);

doaddlabel = false;
if doaddlabel
    labeledge(h, (1 : g.numedges), edgelabels);
end


