%% Pseudo-Two-Dimensional (P2D) Lithium-Ion Battery Model
% This example demonstrates how to setup a P2D model of a Li-ion battery
% and run a simple simulation.

% Clear the workspace and close open figures
clear all
close all

%% Setup the properties of Li-ion battery materials and cell design

jsonstruct = parseBattmoJson('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_lco_graphite_sei.json');

inputparams = BatteryInputParams(jsonstruct);

% We define some shorthand names for simplicity.
elyte   = 'Electrolyte';
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
co      = 'Coating';
am      = 'ActiveMaterial';
am1     = 'ActiveMaterial1';
am2     = 'ActiveMaterial2';
cc      = 'CurrentCollector';
itf     = 'Interface';
sd      = 'SolidDiffusion';
thermal = 'ThermalModel';
ctrl    = 'Control';
sr      = 'SideReaction';
sei     = 'SolidElectrodeInterface';

%% Setup the geometry and computational mesh
% Here, we setup the 1D computational mesh that will be used for the
% simulation. The required discretization parameters are already included
% in the class BatteryGenerator1D. 
gen = BatteryGeneratorP2D();

% Now, we update the inputparams with the properties of the mesh. 
inputparams = gen.updateBatteryInputParams(inputparams);

%%  Initialize the battery model. 
% The battery model is initialized by sending inputparams to the Battery class
% constructor. see :class:`Battery <Battery.Battery>`.
model = GenericBattery(inputparams);

%% The control is used to set up the schedule

jsonstruct.TimeStepping.numberOfTimeSteps = 200;

schedule = model.(ctrl).setupSchedule(jsonstruct);

%% Setup the initial state of the model
% The initial state of the model is setup using the model.setupInitialState() method.

initstate = model.setupInitialState();

%% Setup the properties of the nonlinear solver 
nls = NonLinearSolver(); 
% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10;
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false;
% nls.timeStepSelector=StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);
% Change default tolerance for nonlinear solver
model.nonlinearTolerance = 1e-5*mean(model.Control.ImaxDischarge + model.Control.ImaxCharge);
% Set verbosity
model.verbose = true;

%% Run the simulation
[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

%% Process output and recover the output voltage and current from the output states.
ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
time = cellfun(@(x) x.time, states); 
E    = cellfun(@(x) x.Control.E, states); 
I    = cellfun(@(x) x.Control.I, states);

for istate = 1 : numel(states)
    states{istate} = model.addVariables(states{istate});
end


set(0, 'defaultlinelinewidth', 3)
set(0, 'defaultaxesfontsize', 15)

figure
plot(time/hour, E);
title('Voltage / V')
xlabel('Time / h')

figure
plot(time/hour, I);
title('Current / A')
xlabel('Time / h')

figure
hold on

delta = cellfun(@(state) state.(ne).(co).(am).(sei).delta(end), states);
plot(time/hour, delta/(nano*meter), 'displayname', 'at x_{max}')

delta = cellfun(@(state) state.(ne).(co).(am).(sei).delta(1), states);
plot(time/hour, delta/(nano*meter), 'displayname', 'at x_{min}')

title('SEI thickness in negative electrode/ nm')
xlabel('Time / h')

legend show

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
