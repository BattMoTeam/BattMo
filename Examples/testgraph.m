clear all
close all

% setup mrst modules
mrstModule add ad-core multimodel mrst-gui battery
mrstVerbose off
            
model = Battery_();

model.adminmodel.printPropfunctions;

[g, edgelabels] = setupGraph(model);

plot(g, 'interpreter', 'none');


