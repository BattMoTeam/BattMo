%% Examples of the different control policies
%
% We run the same examples for the different control policies
%
% * |CCDischarge|
% * |CCCharge|
% * |CCCV|
%


%% Setup json input 
%
% We load some parameter sets for the material property and geometry
jsonfilename = fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_lco_graphite.json');

jsonstruct_material = parseBattmoJson(jsonfilename);

jsonfilename = fullfile('Examples', 'JsonDataFiles', 'geometry1d.json');
jsonstruct_geometry = parseBattmoJson(jsonfilename);

%%%
% We merge the material and geometrical parameters in |jsonstruct|

jsonstruct = mergeStructs({jsonstruct_geometry , ...
                               jsonstruct_material
                              });

%% Setup Constant Current Discharge control and run simulation
%

jsonstruct.Control = struct('controlPolicy'     , 'CCDischarge', ...
                            'rampupTime'        , 0.1          , ...
                            'DRate'             , 1            , ...
                            'lowerCutoffVoltage', 2.4);

output = runBattery(jsonstruct);

%%
% We plot the results

set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);
set(0, 'defaultlinelinewidth', 3);

figure
states = output.states;
time = cellfun(@(state) state.time, states);
E    = cellfun(@(state) state.Control.E, states);
plot(time/hour, E);
xlabel('time / h');
ylabel('Voltage / h');


%% Setup Constant Current Charge control and run simulation
%

jsonstruct.Control = struct('controlPolicy'     , 'CCCharge', ...
                            'rampupTime'        , 0.1       , ...
                            'CRate'             , 1         , ...
                            'upperCutoffVoltage', 4.1);

%%%
% We change the state of charge (from 0.99 to 0.01)

jsonstruct.SOC = 0.01;

%%%
% By default the simulation time is set to 1.4/(DRate) hour. We can change that

jsonstruct.TimeStepping.totalTime = 2*hour;

%%%
% We run the simulation

output = runBattery(jsonstruct);


%% Plot Results (1)
%

set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);
set(0, 'defaultlinelinewidth', 3);
figure
states = output.states;
time = cellfun(@(state) state.time, states);
E    = cellfun(@(state) state.Control.E, states);
plot(time/hour, E);
xlabel('time / h');
ylabel('Voltage / h');


%% Setup Constant Current Constant Voltage control and run the simulation
%

jsonstruct.SOC = 0.99;
jsonstruct.Control = struct('controlPolicy'     , 'CCCV'       , ...
                            'initialControl'    , 'discharging', ...
                            'numberOfCycles'    , 2            , ...
                            'CRate'             , 1.5          , ...
                            'DRate'             , 1            , ...
                            'lowerCutoffVoltage', 2.4          , ...
                            'upperCutoffVoltage', 4            , ...
                            'dIdtLimit'         , 2e-6         , ...
                            'dEdtLimit'         , 2e-6);

%%%
% We run the simulation

output = runBattery(jsonstruct);


%% Plot Results (2)

set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);
set(0, 'defaultlinelinewidth', 3);
figure
states = output.states;
time = cellfun(@(state) state.time, states);
E    = cellfun(@(state) state.Control.E, states);
plot(time/hour, E);
xlabel('time / h');
ylabel('Voltage / h');

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
