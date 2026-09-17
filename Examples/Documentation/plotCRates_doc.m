%% Simple plot script used in documentation

close all
figure
hold on

set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);
set(0, 'defaultlinelinewidth', 3);

DRates = [0.8, 1, 2];
for i = 1 : numel(DRates)
    jsonstruct.Control.DRate = DRates(i);
    output = runBattery(jsonstruct);
    plotResult(output);
end

legend({'DRate = 0.8', ...
        'DRate = 1'  , ...
        'DRate = 2'})

axis([0, 1.5, 2.5, 4.5])

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
