set(0, 'defaultlinelinewidth', 3);
set(0, 'defaultaxesfontsize', 15);

% model = ReactionModel();
% model = ReactionModel2();
% model = ThermalModel();
% model = ReactionThermalModel();
% model = UncoupledReactionThermalModel();
% model = ConcentrationModel();
% model = ConcentrationReactionModel();
% model = ConcentrationReactionThermalModel();
model = GenericExampleModel();

cgti = model.cgti;
h = cgti.plot();
% h = cgp.plotModelGraph();
set(h, 'nodefontsize', 20);
% set(h, 'nodefontsize', 9);
set(h, 'linewidth', 5);
set(h, 'arrowsize', 20);

dosave = true;
savedir = '../img';

if dosave
    % filename = 'reacmodelgraph.png';
    % filename = 'reacmodelgraph2.png';
    % filename = 'tempgraph.png';
    % filename = 'concreacgraph.png';
    % filename = 'tempconcreacgraphmodel.png';
    % filename = 'uncoupledreactemp.png';
    filename = 'genericexample.png';
    filename = fullfile(savedir, filename);
    saveas(gcf, filename)
end
    
