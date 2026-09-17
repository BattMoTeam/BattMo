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

cgit = model.cgit;
h = cgit.plot();
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
    

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
