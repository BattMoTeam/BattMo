%% Example of function control
%

%% Case setup
% We use a some standard parameter set for the battery
jsonstruct = parseBattmoJson('Examples/jsondatafiles/sample_input.json');

%% Functional time control
% For a control that depends on time, we switch to a |timeControl| control policy
% 
jsonstruct_control.controlPolicy = 'timeControl';

%%%
% For this control policy, we have to add two functions in the |jsonstruct| structure of the control. The two functions are functions of time. We have
%
% * a |type| function, whose value is either 1, for current, or 2, for voltage.
% * a |value| function, which gives the value of the current or voltage, depending on the type, as a function of time

%% Type control setup
%
% We define the function as a constant equal to 1 for all time. In this case, the value of the control defined below
% will always correspond to a given control current. See |Utilities/JsonSchemas/Function.schema.json| for the syntax of
% a function to be used in json structures.
% 
expression = struct('formula', '1');
jsonstruct_control.type = struct('functionFormat', 'string expression', ...
                                  'argumentList', {{'time'}}, ...
                                  'expression', expression);


%% Value Control setup
%
% We define the function for the control value, in this case a current, because of the type defined above. We use a
% sinusoidal. The |ampere| and |minute| are globally defined variables which can be used here (it would not work with
% local variables).

expression = struct('formula', '1e-2*ampere*sin(2*pi*time/(1*minute))');
jsonstruct_control.value = struct('functionFormat', 'string expression', ...
                                  'argumentList', {{'time'}}, ...
                                  'expression', expression);

%% We finalize the setup of the |jsonstruct| input
%

jsonstruct.Control = jsonstruct_control;

%%%
% We set the total time to 3 minutes

jsonstruct.TimeStepping.totalTime = 3*minute;

%%%
% we set the initial state of charge to 0.5

jsonstruct.SOC = 0.5;

%%%
% we run the simulation

output = runBattery(jsonstruct);

%%
% We plot the results

states = output.states;

time = cellfun(@(state) state.time, states);
E = cellfun(@(state) state.Control.E, states);
I = cellfun(@(state) state.Control.I, states);

figure
yyaxis left
plot(time/minute, E);
ylabel('Voltage / V')
yyaxis right
plot(time/minute, I);
ylabel('Current / A')

%% Tabulated data for the function control
%
% We change to tabulated data for the current control value 

dataX = [0, 1*minute, 2*minute, 3*minute];
dataY = ampere*[0, 1e-2, 1e-2, 0];
expression = struct('formula', '1e-2*ampere*sin(2*pi*time/(1*minute))');
jsonstruct_control.value = struct('functionFormat', 'tabulated', ...
                                  'argumentList'  , {{'time'}} , ...
                                  'dataX'         , dataX      , ...
                                  'dataY'         , dataY);

%%
% We re-run the simulation

jsonstruct.Control = jsonstruct_control;

output = runBattery(jsonstruct);

%%%
% We plot the results
%

states = output.states;

time = cellfun(@(state) state.time, states);
E = cellfun(@(state) state.Control.E, states);
I = cellfun(@(state) state.Control.I, states);

figure
yyaxis left
plot(time/minute, E);
ylabel('Voltage / V')
yyaxis right
plot(time/minute, I);
ylabel('Current / A')

%% Tabulated current time series
%
% A measured or otherwise prescribed current time series can be supplied
% directly as tabulated data. The first column is time and the second is
% current. The current is linearly interpolated between the sample points.
% Positive current corresponds to discharge and negative current to charge.

currentTimeSeries = [0.00,  0.0; ...
                     0.25,  1.0; ...
                     1.00,  1.0; ...
                     1.25,  0.0; ...
                     2.00, -0.5; ...
                     2.75, -0.5; ...
                     3.00,  0.0];

currentTimeSeries(:, 1) = currentTimeSeries(:, 1)*minute;
currentTimeSeries(:, 2) = currentTimeSeries(:, 2)*1e-2*ampere;

jsonstruct_control.type = struct('functionFormat', 'constant', ...
                                 'argumentList'  , {{'time'}}, ...
                                 'value'         , 1);

jsonstruct_control.value = struct('functionFormat', 'tabulated', ...
                                  'argumentList'  , {{'time'}}, ...
                                  'dataX'         , currentTimeSeries(:, 1), ...
                                  'dataY'         , currentTimeSeries(:, 2));

jsonstruct.Control                = jsonstruct_control;
jsonstruct.TimeStepping.totalTime = currentTimeSeries(end, 1);

output = runBattery(jsonstruct);

states = output.states;

time = cellfun(@(state) state.time, states);
E = cellfun(@(state) state.Control.E, states);
I = cellfun(@(state) state.Control.I, states);

figure
yyaxis left
plot(time/minute, E, 'DisplayName', 'Voltage');
ylabel('Voltage / V')
yyaxis right
plot(time/minute, I, 'DisplayName', 'Simulated current');
hold on
plot(currentTimeSeries(:, 1)/minute, ...
     currentTimeSeries(:, 2)/ampere, 'o', ...
     'DisplayName', 'Current samples');
hold off
ylabel('Current / A')
xlabel('Time / min')
legend('Location', 'best')

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
