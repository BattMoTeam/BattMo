%% Example with non-linear solver

%% Setup nonlinear solver

jsonfilename = '/home/xavier/Matlab/Projects/battmo/Examples/JsonDataFiles/linear_solver_setup.json';
jsonstruct_nonlinear_solver = parseBattmoJson(jsonfilename);

%% Setup material properties

jsonfilename = fullfile('ParameterData'        , ...
                        'BatteryCellParameters', ...
                        'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_nmc_graphite.json');
jsonstruct_material = parseBattmoJson(jsonfilename);

jsonstruct_material = removeJsonStructFields(jsonstruct_material           , ...
                                             {'include_current_collectors'}, ...
                                             {'ThermalModel', 'externalHeatTransferCoefficient'});

%% Setup geometry
jsonfilename = fullfile('Examples'     , ...
                        'JsonDataFiles', ...
                        'geometry3d.json');
jsonstruct_geometry = parseBattmoJson(jsonfilename);

%% Setup Control
jsonfilename = fullfile('Examples', 'JsonDataFiles', 'cc_discharge_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);


jsonstruct = mergeJsonStructs({jsonstruct_nonlinear_solver, ...
                               jsonstruct_geometry        , ...
                               jsonstruct_material        , ...
                               jsonstruct_control});


%% Run the simulation

output = runBatteryJson(jsonstruct);

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

Tmin = cellfun(@(state) min(state.ThermalModel.T + T0), states);
Tmax = cellfun(@(state) max(state.ThermalModel.T + T0), states);

figure
hold on
plot(time / hour, Tmin, 'displayname', 'min T');
plot(time / hour, Tmax, 'displayname', 'max T');
title('Temperature / C')
xlabel('time / h');
ylabel('Temperature / C');

legend show

%% change setup for cooling coefficient

jsonstruct.ThermalModel.externalHeatTransferCoefficientTab = 100;
jsonstruct.ThermalModel.externalHeatTransferCoefficient = 0;

output = runBatteryJson(jsonstruct);

%% Plot the minimum and maximum values of the temperature

states = output.states;

Tmin = cellfun(@(state) min(state.ThermalModel.T + T0), states);
Tmax = cellfun(@(state) max(state.ThermalModel.T + T0), states);
time = output.time;

figure
hold on
plot(time / hour, Tmin, 'displayname', 'min T');
plot(time / hour, Tmax, 'displayname', 'max T');
title('Temperature / C')
xlabel('time / h');
ylabel('Temperature / C');

legend show

%% plot final temperature distribution

state = states{end}
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
