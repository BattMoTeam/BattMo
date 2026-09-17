clear
close all
clc

jsonstruct = parseBattmoJson('params_initial.json');

% Control
jsonfilename = fullfile('Examples', 'JsonDataFiles', 'cccv_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);

jsonstruct_control.Control.lowerCutoffVoltage = 4;
jsonstruct_control.Control.upperCutoffVoltage = 4.8;
jsonstruct_control.Control.DRate              = 1;
jsonstruct_control.Control.CRate              = 1;
jsonstruct_control.Control.dIdtLimit          = 1E-4;
jsonstruct_control.Control.dEdtLimit          = 1E-2;
jsonstruct_control.Control.numberOfCycles     = 5;

jsonstruct = mergeStructs({jsonstruct_control, jsonstruct});

jsonstruct.TimeStepping.numberOfTimeSteps = 1000;

% For this simulation, to resolve the control switch we need to
% increase the number of time step cut allowed.
jsonstruct.NonLinearSolver.maxTimestepCuts = 20;

output = runBattery(jsonstruct, 'runSimulation', true);

% Integrate the current over time to calculate capacity
capacity = cumtrapz(output.time, output.I);

% Ensure the capacity starts at 0 and reverses back to 0
% capacity = capacity - min(capacity); % Shift to start from 0
% capacity = max(capacity) - capacity; % Reverse back to 0 when fully charged

figure;
plot(capacity / hour, output.E, 'linewidth', 2);
xlabel('Capacity  /  A \cdot h')
ylabel('Cell Voltage  /  V')
set(gca, 'FontSize', 18);

figure
plot(output.time / hour, output.E, 'linewidth', 2);
xlabel('Time  /  h')
ylabel('Cell Voltage  /  V')
set(gca, 'FontSize', 18);

plotContours(output)

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
