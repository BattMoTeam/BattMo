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

writeJsonStruct(jsonstruct, 'test.json');

jsonstruct_control = jsonstruct;

jsonstruct_material = parseBattmoJson(fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', 'lithium_ion_battery_nmc_graphite.json'));
jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));

jsonstruct_material = removeJsonStructFields(jsonstruct_material              , ...
                                             {'Control', 'DRate'}             , ...
                                             {'Control', 'controlPolicy'}     , ...
                                             {'Control', 'upperCutoffVoltage'}, ...
                                             {'Control', 'rampupTime'}        , ...
                                             {'Control', 'lowerCutoffVoltage'});

jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                               jsonstruct_geometry, ...
                               jsonstruct_control});

output = runBatteryJson(jsonstruct);

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
