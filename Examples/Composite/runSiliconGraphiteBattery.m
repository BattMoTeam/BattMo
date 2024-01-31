%% Composite Silicon Graphite electrode
clear
close all

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core

%% Shortcuts
% We define shorcuts for the sub-models.

ne   = 'NegativeElectrode';
pe   = 'PositiveElectrode';
co   = 'Coating';
am1  = 'ActiveMaterial1';
am2  = 'ActiveMaterial2';
bd   = 'Binder';
ad   = 'ConductingAdditive';
sd   = 'SolidDiffusion';
itf  = 'Interface';
ctrl = 'Control';

%% Setup the properties of the battery
%
% We load the property of a composite silicon graphite electrode, see
% :ref:`compositeElectrode`
%

jsonstruct_composite_material = parseBattmoJson('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_silicon_graphite.json');

%%
% For the remaining properties, we consider a standard data set
jsonstruct_cell = parseBattmoJson('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json');

%%
% We remove form the standard data set the :code:`ActiveMaterial`
% field. This step is not necessary but is cleaner and we avoid a
% warning.
jsonstruct_cell.(ne).(co) = rmfield(jsonstruct_cell.(ne).(co), 'ActiveMaterial');

%%
% We merge the two json structures
jsonstruct = mergeJsonStructs({jsonstruct_composite_material, ...
                               jsonstruct_cell});

%%
% We do not consider the thermal model and remove the current
% collector
jsonstruct.use_thermal = false;
jsonstruct.include_current_collectors = false;

%%
% We instantiate the battery :code:`InputParams` object
inputparams = BatteryInputParams(jsonstruct);

%%
% We set the mass fractions of the different material in the coating
% of the negative electrode. This information could have been passed
% in the json file earlier (:ref:`compositeElectrode`)

inputparams.(ne).(co).(am1).massFraction = 0.9;
inputparams.(ne).(co).(am2).massFraction = 0.08;
inputparams.(ne).(co).(bd).massFraction  = 0.01;
inputparams.(ne).(co).(ad).massFraction  = 0.01;


%%
% Now, we update the inputparams with the properties of the grid.
gen = BatteryGeneratorP2D();
inputparams = gen.updateBatteryInputParams(inputparams);

%% We change the given CRate

CRate = 0.1;
inputparams.(ctrl).CRate = CRate;

% Allow for switching to voltage control
inputparams.(ctrl).useCVswitch = true;

%% Model Instantiation
% We instantiate the model

model = Battery(inputparams);


%% Setup schedule (control and time stepping)
% We will simulate two consecutive periods: a discharge followed by a
% charge.
%
% We start with the charge period

step    = model.(ctrl).setupScheduleStep();
control = model.(ctrl).setupScheduleControl();

schedule = struct('control', control, 'step', step);

%% Setup the initial state of the model
%
% We use the default initialisation given by a method in the model

initstate = model.setupInitialState();

%% Setup the properties of the nonlinear solver
% We adjust some settings for the nonlinear solver
nls = NonLinearSolver();

%%
% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10;
%%
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false;
%%
% We use a time step selector based on absolute change of a target
% value, in our case the output voltage
nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, ...
                                                   'targetChangeAbs', 0.1);
%%
% We adjust the nonlinear tolerance
model.nonlinearTolerance = 1e-3*model.Control.Imax;

%% Run the simulation for the discharge

fprintf('Run the simulation for the discharge\n');
[~, states] = simulateScheduleAD(initstate, model, schedule, ...
                                 'OutputMinisteps', true, ...
                                 'NonLinearSolver', nls);

dischargeStates = states;

%% Setup charge schedule

% We use the last computed state of the discharge as the initial state
% for the charge period.
initstate = states{end};

% We use a new control model.

inputparams.Control = CCChargeControlModelInputParams(jsonstruct.Control);
inputparams.(ctrl).CRate = CRate;

inputparams = inputparams.validateInputParams();

model = Battery(inputparams);

step    = model.(ctrl).setupScheduleStep();
control = model.(ctrl).setupScheduleControl();

% Use the control in the schedule
schedule = struct('control', control, 'step', step);

%% Run the simulation for the charge period
fprintf('Run the simulation for the charge\n');
[~, states] = simulateScheduleAD(initstate, model, schedule, ...
                                 'OutputMinisteps', true, ...
                                 'NonLinearSolver', nls);

chargeStates = states;

%% Visualisation

%%
% We concatenate the states we have computed
allStates = vertcat(dischargeStates, chargeStates);

%%
% We extract the voltage, current and time from the simulation output
ind = cellfun(@(x) ~isempty(x), allStates);
allStates = allStates(ind);
E    = cellfun(@(x) x.Control.E, allStates);
I    = cellfun(@(x) x.Control.I, allStates);
time = cellfun(@(x) x.time, allStates);

%%
%  We plot the voltage and current
figure
subplot(2, 1, 1);
plot(time/hour, E);
xlabel('Time / h');
ylabel('Voltage / V');
title('Voltage')
subplot(2, 1, 2);
plot(time/hour, I/milli);
xlabel('Time / h');
ylabel('Current / mA');
title('Current')

%%
% We compute and plot the state of charges in the different material

figure
hold on

for istate = 1 : numel(allStates)
    allStates{istate} = model.evalVarName(allStates{istate}, {ne, co, 'SOC'});
end

SOC  = cellfun(@(x) x.(ne).(co).SOC, allStates);
SOC1 = cellfun(@(x) x.(ne).(co).(am1).SOC, allStates);
SOC2 = cellfun(@(x) x.(ne).(co).(am2).SOC, allStates);

plot(time/hour, SOC, 'displayname', 'SOC - cumulated');
plot(time/hour, SOC1, 'displayname', 'SOC - Graphite');
plot(time/hour, SOC2, 'displayname', 'SOC - Silicon');

xlabel('Time / h');
ylabel('SOC / -');
title('SOCs')

legend show



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
