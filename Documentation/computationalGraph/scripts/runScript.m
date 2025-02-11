set(0, 'defaultlinelinewidth', 3);
set(0, 'defaultaxesfontsize', 15);

% model = ReactionModel();
% model = ReactionModel2();
% model = ThermalModel();
% model = ReactionThermalModel();
model = UncoupledReactionThermalModel();
% model = ConcentrationModel();
% model = ConcentrationReactionModel();
% model = ConcentrationReactionThermalModel();
cgp = model.cgp;
cgt = model.cgt;
h = cgp.plot();
% h = cgp.plotModelGraph();
set(h, 'nodefontsize', 14);
% set(h, 'nodefontsize', 9);
set(h, 'linewidth', 3);
set(h, 'arrowsize', 20);

dosave = true;
savedir = '../img';

if dosave
    % filename = 'reacmodelgraph.png';
    % filename = 'reacmodelgraph2.png';
    % filename = 'tempgraph.png';
    % filename = 'concreacgraph.png';
    % filename = 'tempconcreacgraphmodel.png';
    filename = 'uncoupledreactemp.png';
    filename = fullfile(savedir, filename);
    saveas(gcf, filename)
end
    
