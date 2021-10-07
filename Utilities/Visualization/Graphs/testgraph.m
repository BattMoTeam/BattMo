clear all
close all

% setup mrst modules
mrstModule add ad-core multimodel mrst-gui battery 
mrstVerbose off

set(0, 'DefaultAxesFontSize', 16);
set(0, 'defaulttextfontsize', 18);
set(0, 'DefaultFigurePosition', [278 109 1570 822]);

gstyle = {'interpreter', 'none', 'LineWidth', 3, 'ArrowSize', 20, 'NodeFontSize', 14};

%% Choose the model you want to use below by uncommenting it and commenting the others.

% model = ActiveMaterial_('am');
% model = ElectrodeActiveComponent_('eac');
% model = ElectronicComponent_('eac');
% model = ElectroChemicalComponent_('eac');
% model = CurrentCollector_('cc');
model = Electrolyte_('elyte');
% model = orgLiPF6_('cc');
% model = Electrode_('pn');
% model = Battery_();

% setup the model
model = model.initiateCompositeModel();

% plot the computational graph
[g, edgelabels] = setupGraph(model);
plot(g, gstyle{:});

% print the functions with input and output arguments
model.adminmodel.printPropfunctions;

% print the functions in call order (such as all the computational graph is spanned)
% This require external routines that can be found in the MATLAB BGL code base (see https://se.mathworks.com/matlabcentral/fileexchange/10922-matlabbgl)
doorder = false;
try 
    mrstModule add bgl
    doorder = true;
catch
    fprintf('To order the graph function you need to install matlab bgl module.\n see https://se.mathworks.com/matlabcentral/fileexchange/10922-matlabbgl');
end
if doorder
    fprintf('\n\n\n**** Functions ordered by argument dependancy:\n\n');
    setupOrderedFunctions(model);
end







