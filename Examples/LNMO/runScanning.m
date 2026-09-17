clear all

sc = Scanning();

thickness = linspace(50*micro*meter, 100*micro*meter, 5);
SiContent = linspace(0.01, 0.1, 5);

[inds, inputs] = cartesianProduct(thickness, SiContent);

for iinput = 1 : numel(inputs)
    
    input = inputs{iinput};
    sc.computeSpecificEnergy(input.thickness, input.SiContent);

    E(iinput) = sc.specificEnergyHour;

end


%%

close all

set(0, 'defaultlinelinewidth', 3);
set(0, 'defaultaxesfontsize', 15);

X = repmat(thickness', 1, numel(SiContent))/(micro*meter);
Y = repmat(SiContent, numel(thickness), 1);
Z = reshape(E, numel(thickness), numel(SiContent));

% figure
% surf(X, Y, Z);
% xlabel('Negative Electrode Length / μm');
% ylabel('Silicon Ratio (mass fraction)');
% zlabel('Specific Energy / Wh \cdot kg^{-1}');

figure

contourf(X, Y, Z);
xlabel('Negative Electrode Length / μm');
ylabel('Silicon Ratio (mass fraction)');
zlabel('Specific Energy / Wh \cdot kg^{-1}');

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
