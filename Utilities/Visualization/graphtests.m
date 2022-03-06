
clear all
close all

%% Add MRST module
mrstModule add ad-core

modelcase = 'battery';

switch modelcase
  case 'elyte'
    G = cartGrid([10, 10]);
    G = computeGeometry(G);
    model = orgLiPF6('elyte', G, (1 : G.cells.num));
  case 'battery'
    model = Battery();
    model.I = 0.1;
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





%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BatMo

BatMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BatMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BatMo.  If not, see <http://www.gnu.org/licenses/>.
%}
