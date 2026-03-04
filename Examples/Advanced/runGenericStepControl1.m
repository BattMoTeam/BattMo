%% Example using generic control input
%

%% input setup
%

%%
% We load material parameter input
jsonstruct_material = parseBattmoJson(fullfile('ParameterData', ...
                                               'BatteryCellParameters', ...
                                               'LithiumIonBatteryCell', ...
                                               'lithium_ion_battery_nmc_graphite.json'));

%%
% We load geometry parameter input


jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));

%%
% We load the generic control parameter we want to use

jsonstruct_control  = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'generic_step_control_example1.json'));

%%
% We print it to command window

viewJsonStruct(jsonstruct_control)

%%
% We remove the fields that are not relevant for this example. This is done to simplify the input parameters and focus
% on the relevant aspects of the model.
%

jsonstruct_material = removeJsonStructFields(jsonstruct_material              , ...
                                             {'Control', 'DRate'}             , ...
                                             {'Control', 'controlPolicy'}     , ...
                                             {'Control', 'upperCutoffVoltage'}, ...
                                             {'Control', 'rampupTime'}        , ...
                                             {'Control', 'lowerCutoffVoltage'});

%%
% We merge all the input structures
% 
jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                               jsonstruct_geometry, ...
                               jsonstruct_control});


%%
% We do not include thermal effects and the current collectors in this example to simplify the model.
%
jsonstruct.use_thermal = false;
jsonstruct.include_current_collectors = false;

output = runBatteryJson(jsonstruct);

%% Process output and recover the output voltage and current from the output states.

states = output.states;

states = states(ind);
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
