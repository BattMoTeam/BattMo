clear all
close all

%% Example with non-linear solver

%% Setup nonlinear solver

jsonFilename = fullfile(battmoDir(), 'Examples', 'JsonDataFiles', 'linear_solver_setup.json');
solverJson = parseBattmoJson(jsonFilename);

%% Setup material properties

jsonFilename = fullfile(battmoDir(), 'ParameterData', 'BatteryCellParameters', ...
    'LithiumIonBatteryCell', 'lithium_ion_battery_lco_graphite.json');
materialJson = parseBattmoJson(jsonFilename);

% Use the current collectors and cooling coefficients specified by the geometry.
materialJson = removeStructFields(materialJson, ...
    {'include_current_collectors'}, ...
    {'ThermalModel', 'externalHeatTransferCoefficient'});

%% Setup geometry

jsonFilename = fullfile(battmoDir(), 'Examples', 'JsonDataFiles', 'geometry3d.json');
geometryJson = parseBattmoJson(jsonFilename);

%% Setup control

jsonFilename = fullfile(battmoDir(), 'Examples', 'JsonDataFiles', 'cc_discharge_control.json');
controlJson = parseBattmoJson(jsonFilename);

jsonstruct = mergeStructs({solverJson, geometryJson, materialJson, controlJson});

%% Run the simulation

output = runBattery(jsonstruct);

%% plot voltage

time = output.time;
E    = output.E;

set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);
set(0, 'defaultlinelinewidth', 3);

figure
plot(time/hour, E)
title('Voltage / V');
xlabel('time / h');
ylabel('voltage / V');

%% Plot the minimum and maximum values of the temperature

T0 = PhysicalConstants.absoluteTemperature;

states = output.states;

minTemperature = cellfun(@(state) min(state.ThermalModel.T + T0), states);
maxTemperature = cellfun(@(state) max(state.ThermalModel.T + T0), states);

figure
hold on
plot(time / hour, minTemperature, 'displayname', 'min T');
plot(time / hour, maxTemperature, 'displayname', 'max T');
title('Temperature / C')
xlabel('time / h');
ylabel('Temperature / C');

legend show

%% change setup for cooling coefficient

jsonstruct.ThermalModel.externalHeatTransferCoefficientTab = 100;
jsonstruct.ThermalModel.externalHeatTransferCoefficient = 0;

output = runBattery(jsonstruct);

%% Plot the minimum and maximum values of the temperature

states = output.states;

minTemperature = cellfun(@(state) min(state.ThermalModel.T + T0), states);
maxTemperature = cellfun(@(state) max(state.ThermalModel.T + T0), states);
time = output.time;

figure
hold on
plot(time / hour, minTemperature, 'displayname', 'min T');
plot(time / hour, maxTemperature, 'displayname', 'max T');
title('Temperature / C')
xlabel('time / h');
ylabel('Temperature / C');

legend show

%% plot final temperature distribution

model = output.model;
state = states{end};
figure
plotCellData(model.ThermalModel.grid, ...
    state.ThermalModel.T + T0);
colorbar
title('Temperature / C');
view([50, 50]);

%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
