clear all
close all

printer = @(s) disp(jsonencode(s, 'PrettyPrint', true));

experiment = {
    'Rest for 4000 s';
    'Discharge at 500 mA until 3.0 V';
    'Hold at 3.0 V until 1e-4 A';
    'Charge at 1 A until 4.0 V';
    'Rest for 1 hour';
             };

jsonstruct = convertPyBaMMtoJson(experiment);

printer(jsonstruct);

writeStruct(jsonstruct, 'test.json');

jsonstruct_control = jsonstruct;

jsonstruct_material = parseBattmoJson(fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', 'lithium_ion_battery_lco_graphite.json'));
jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));

jsonstruct_material = removeStructFields(jsonstruct_material              , ...
                                             {'Control', 'DRate'}             , ...
                                             {'Control', 'controlPolicy'}     , ...
                                             {'Control', 'upperCutoffVoltage'}, ...
                                             {'Control', 'rampupTime'}        , ...
                                             {'Control', 'lowerCutoffVoltage'});

jsonstruct = mergeStructs({jsonstruct_material, ...
                               jsonstruct_geometry, ...
                               jsonstruct_control});

output = runBattery(jsonstruct);

%% Process output and recover the output voltage and current from the output states.

states = output.states;

E = cellfun(@(x) x.Control.E, states);
I = cellfun(@(x) x.Control.I, states);
time = cellfun(@(x) x.time, states);

figure
plot(time/hour, E, '*-', 'linewidth', 3);
grid on
xlabel 'time  / h';
ylabel 'potential  / V';

figure
plot(time/hour, I, '*-', 'linewidth', 3);
grid on
xlabel 'time  / h';
ylabel 'Current  / A';

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
