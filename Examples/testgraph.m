clear all
close all

% setup mrst modules
mrstModule add ad-core multimodel mrst-gui battery
mrstVerbose off
            
model = Battery_();
% model = CurrentCollector_('cc');
% model = model.initiateCompositeModel();

model.adminmodel.printPropfunctions;

[g, edgelabels] = setupGraph(model);

plot(g, 'interpreter', 'none');

f = fopen('test.csv', 'w');
edges = g.Edges;
edges = edges.Variables;
for ind = 1 : size(edges, 1)
    fprintf(f, '%s;%s\n', edges{ind, 1}, edges{ind, 2});
end
fclose(f);



