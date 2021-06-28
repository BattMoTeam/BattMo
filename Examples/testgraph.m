clear all
close all

% setup mrst modules
mrstModule add ad-core multimodel mrst-gui battery
mrstVerbose off

set(0, 'DefaultAxesFontSize', 16);
set(0, 'defaulttextfontsize', 18);
set(0, 'DefaultFigurePosition', [278 109 1570 822]);

gstyle = {'interpreter', 'none', 'LineWidth', 3, 'ArrowSize', 20, 'NodeFontSize', 14};

% model = ActiveMaterial_('am');
% model = ElectrodeActiveComponent_('eac');
% model = ElectronicComponent_('eac');
% model = ElectroChemicalComponent_('eac');
% model = CurrentCollector_('cc');
model = Electrolyte_('cc');
% model = Electrode_('pn');
% model = Electrolyte_('elyte');
% model = Battery_();

model = model.initiateCompositeModel();
model.adminmodel.printPropfunctions;

[g, edgelabels] = setupGraph(model);

plot(g, gstyle{:});

f = fopen('test.sif', 'w');
edges = g.Edges;
edges = edges.Variables;
for ind = 1 : size(edges, 1)
    fprintf(f, '%s x %s\n', edges{ind, 1}, edges{ind, 2});
end
fclose(f);



