clear all

jsonstruct = parseBattmoJson('Examples/jsondatafiles/sample_input.json');

jsonstruct_control.controlPolicy = 'timeControl';

expression = struct('formula', '1e-2*ampere*sin(2*pi*time/(1*minute))');
jsonstruct_control.value = struct('functionFormat', 'string expression', ...
                                  'argumentList', {{'time'}}, ...
                                  'expression', expression);


dataX = [0, 1*minute, 2*minute, 3*minute];
dataY = ampere*[0, 1e-2, 1e-2, 0];
expression = struct('formula', '1e-2*ampere*sin(2*pi*time/(1*minute))');
jsonstruct_control.value = struct('functionFormat', 'tabulated', ...
                                  'argumentList'  , {{'time'}} , ...
                                  'dataX'         , dataX      , ...
                                  'dataY'         , dataY);

% type control is either 1 (current) or 2 (voltage), see Control/TimeControlModel.m
expression = struct('formula', '1');
jsonstruct_control.type = struct('functionFormat', 'string expression', ...
                                  'argumentList', {{'time'}}, ...
                                  'expression', expression);


jsonstruct.Control = jsonstruct_control;

jsonstruct.TimeStepping.totalTime = 3*minute;

jsonstruct.SOC = 0.5;

output = runBatteryJson(jsonstruct, 'runSimulation', true);

%%

states = output.states;

time = cellfun(@(state) state.time, states);
E = cellfun(@(state) state.Control.E, states);
I = cellfun(@(state) state.Control.I, states);

close all

figure
yyaxis left
plot(time/minute, E);
ylabel('Voltage / V')
yyaxis right
plot(time/minute, I);
ylabel('Current / A')
