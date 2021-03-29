clear all
close all

% setup mrst modules
mrstModule add ad-core multimodel mrst-gui battery matlab_bgl
mrstVerbose off
      
set(0, 'DefaultAxesFontSize', 16);
set(0, 'defaulttextfontsize', 18);
set(0, 'DefaultFigurePosition', [278 109 1570 822]);

%% ElectronicComponent

gstyle = {'interpreter', 'none', 'LineWidth', 3, 'ArrowSize', 20, 'NodeFontSize', 18};

figure
model = ElectronicComponent_('main');
model = model.initiateCompositeModel();
[g, edgelabels] = setupGraph(model);
plot(g, gstyle{:})
title('ElectronicComponent');

model.adminmodel.printPropfunctions;


%% ElectroChemicalComponent

gstyle = {'interpreter', 'none', 'LineWidth', 3, 'ArrowSize', 20, 'NodeFontSize', 14};

figure
model = ElectroChemicalComponent_('main');
model = model.initiateCompositeModel();
[g, edgelabels] = setupGraph(model);
plot(g, gstyle{:});
title('ElectroChemicalComponent');

model.adminmodel.printPropfunctions;


%% ActiveMaterial

figure
model = ActiveMaterial_('main');
model = model.initiateCompositeModel();
[g, edgelabels] = setupGraph(model);
plot(g, gstyle{:});
title('ActiveMaterial');

model.adminmodel.printPropfunctions;

%% ActiveElectroChemicalComponent

figure
model = ActiveElectroChemicalComponent_('main');
model = model.initiateCompositeModel();
[g, edgelabels] = setupGraph(model);
plot(g, gstyle{:});
title('ActiveElectroChemicalComponent');

model.adminmodel.printPropfunctions;


%% Electrode

figure
model = Electrode_('main');
model = model.initiateCompositeModel();
[g, edgelabels] = setupGraph(model);
plot(g, gstyle{:});
title('Electrode');

model.adminmodel.printPropfunctions;

%% Battery

figure
model = Battery_();
model = model.initiateCompositeModel();
[g, edgelabels] = setupGraph(model);
plot(g, gstyle{:});
title('Battery');

model.adminmodel.printPropfunctions;


setupOrderedGraph(model, 'temp.txt');

return

dosave = false;
if dosave
    f = fopen('test.csv', 'w');
    edges = g.Edges;
    edges = edges.Variables;
    for ind = 1 : size(edges, 1)
        fprintf(f, '%s;%s\n', edges{ind, 1}, edges{ind, 2});
    end
    fclose(f);
end


