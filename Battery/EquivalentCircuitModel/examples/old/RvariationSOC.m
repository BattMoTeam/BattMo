%% Equivalent Circuit Model Simulation

% Load data set

params = createParametersECM();
inputparams = EquivalentCircuitModelInputParams(params);

% Setup model
model = EquivalentCircuitModel(inputparams);

% Run Simulation
[t, U, I, SOC] = model.solve();


% Plotting
figure(Name='Highlighting the difference in resistance depending on the state of charge'); 
tiledlayout(3, 1);

nexttile
plot(t, U, 'LineWidth', 2)
title('Battery Voltage')
xlabel('Time /s')
ylabel('Voltage /V')
grid on

nexttile
plot(t, I, 'LineWidth', 2, 'Color', 'r')
title('Applied Current')
xlabel('Time /s')
ylabel('Current /A')
grid on

nexttile
plot(t, SOC, 'LineWidth', 2, 'Color', 'g')
title('SOC')
xlabel('Time /s')
ylabel('SOC')
grid on

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
