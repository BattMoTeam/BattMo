clear all

dt = 1*second;
T  = 3*minute;
N  = T/dt;

step.val = dt*ones(N, 1);
step.control = ones(N, 1);

period = 1*minute;
control.src = @(time) (1e-2*ampere*sin(2*pi*time/period));
control.CC = true;

schedule.step = step;
schedule.control = control;

jsonstruct = parseBattmoJson('Examples/jsondatafiles/sample_input.json');

jsonstruct.Control = [];
jsonstruct.Control.controlPolicy = 'CC';

jsonstruct.SOC = 0.5;

model = setupModelFromJson(jsonstruct);

initstate = model.setupInitialState();

[~, states] = simulateScheduleAD(initstate, model, schedule);

%%

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
