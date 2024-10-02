%% Example: Galvanostatic Intermittent Charge Titration Technique (GITT) Simulation of a Lithium-Ion Battery Model
% This example demonstrates how to apply a GITT protocol in a Li-ion
% battery simulation.

% Clear the workspace and close open figures
clear
close all


%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa
mrstVerbose off

% We define some shorthand names for simplicity.
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
cc      = 'CurrentCollector';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
ctrl    = 'Control';

%% Setup the properties of Li-ion battery materials and cell design
% The properties and parameters of the battery cell, including the
% architecture and materials, are set using an instance of
% :class:`BatteryInputParams <Battery.BatteryInputParams>`. This class is
% used to initialize the simulation and it propagates all the parameters
% throughout the submodels. The input parameters can be set manually or
% provided in json format. All the parameters for the model are stored in
% the inputparams object.
jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));

jsonstruct.include_current_collectors = true;

inputparams = BatteryInputParams(jsonstruct);


%% Setup the geometry and computational grid
% Here, we setup the geometry and computational grid that will be used for
% the simulation. The user can select the dimensionality of the model. The
% required discretization parameters are already included in the associated
% battery generator class.
modelcase = '2D';

% Generate the battery based on the selected dimensionality in modelcase
switch modelcase
  case '1D'

    gen = BatteryGeneratorP2D();

  case '2D'
    gen = BatteryGeneratorP3D();

    % To avoid convergence issues, override the electronic conductivity
    inputparams.(ne).(cc).effectiveElectronicConductivity = 0.1*inputparams.(ne).(cc).electronicConductivity;
    inputparams.(pe).(cc).effectiveElectronicConductivity = 0.1*inputparams.(pe).(cc).electronicConductivity;

  case '3D'
    gen = BatteryGeneratorP4D();

    fac = 1;
    gen.facx = fac;
    gen.facy = fac;
    gen.facz = fac;
    gen = gen.applyResolutionFactors();

end

inputparams = gen.updateBatteryInputParams(inputparams);

inputparams.(ctrl).useCVswitch = true;


%%  Initialize the battery model.
% The battery model is initialized by sending inputparams to the Battery class
% constructor. see :class:`Battery <Battery.Battery>`.
model = GenericBattery(inputparams);

%% Plot
plotBatteryGrid(model, 'setstyle', false);

%% Compute the nominal cell capacity and choose a C-Rate
% The nominal capacity of the cell is calculated from the active materials.
% This value is then combined with the user-defined C-Rate to set the cell
% operational current.
C      = computeCellCapacity(model);
DRate  = 2;
inputI = (C/hour)*DRate;
inputE = jsonstruct.Control.upperCutoffVoltage;

%% Setup the parameters of the GITT protocol
pulseFraction  = 0.01;
relaxationTime = 4*hour;
switchTime     = 1*milli*second; % switching time (linear interpolation between the two states)
dischargeTime  = pulseFraction*DRate*hour; % time of discharging

% Discretization parameters
tfac                          = 1;
intervalsPerGalvanostaticStep = 5*tfac; % Number of time step in galvanostatic phase
intervalsPerRelaxationStep    = 5*tfac; % Number of time step in relaxation phase
intervalsPerRampupStep        = 3*tfac; % Number of time step in rampup phase

testing = true;
if testing
    fprintf('We setup a smaller case for quicker testing\n');
    numberOfIntervals = 3;
end


%% Setup the initial state of the model
% The initial state of the model is dispatched using the
% model.setupInitialState() method.
initstate = model.setupInitialState();

%% Setup the time step schedule
% Smaller time steps are used to ramp up the current from zero to its
% operational value. Larger time steps are then used for the normal
% operation.
time_init = 5*minute;
dtpoints = [switchTime; ...
            dischargeTime - switchTime;
            switchTime;
            relaxationTime - switchTime];
dIpoints = [1; 0; -1; 0];

tpoints = [0; cumsum(repmat(dtpoints, numberOfIntervals, 1))];
Ipoints = [0; cumsum(repmat(dIpoints, numberOfIntervals, 1))];

tpoints = time_init + tpoints; % starts at end of activation phase
Ipoints = inputI*Ipoints; % scales with Iinput

%% Setup the operating limits for the cell
% The maximum and minimum voltage limits for the cell are defined using
% stopping and source functions. A stopping function is used to set the
% lower voltage cutoff limit. A source function is used to set the upper
% voltage cutoff limit.

tup = 1*milli*second; % rampup time
srcfunc_init = @(time, I, E, inputI) rampupSwitchControl(time, tup, I, E, inputI, inputE);
srcfunc_gitt = @(time, I, E, dummy) tabulatedIControl(time, tpoints, Ipoints);

control = model.Control.setupScheduleControl();

control(1)     = control;
control(1).src = srcfunc_init;
control(2)     = control;
control(2).src = srcfunc_gitt;

n_init = 5;
dt_init = rampupTimesteps(time_init, time_init/n_init, 3);

dt_cycle = [switchTime/intervalsPerRampupStep*ones(intervalsPerRampupStep, 1); ...
            (dischargeTime - switchTime)/intervalsPerGalvanostaticStep*ones(intervalsPerGalvanostaticStep, 1); ...
            switchTime/intervalsPerRampupStep*ones(intervalsPerRampupStep, 1); ...
            (relaxationTime - switchTime)/intervalsPerRelaxationStep*ones(intervalsPerRelaxationStep, 1)];

dt_cycle = repmat(dt_cycle, numberOfIntervals, 1);

step.val = [dt_init; dt_cycle];
step.control = [ones(numel(dt_init), 1); ...
                2*ones(numel(dt_cycle), 1)];

schedule = struct('control', control, 'step', step);

%% Setup nonlinear solver
nls = NonLinearSolver();

% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10;
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false;
% Change default tolerance for nonlinear solver
model.nonlinearTolerance = 1e-4;
% Get more or less verbose output
model.verbose = true;

use_amg = false;
if use_amg
    mrstModule add agmg
    nls.LinearSolver = LinearSolverBattery('method', 'agmg', 'verbosity', 0);
end

%% Run simulation
doprofiling = false;
if doprofiling
    profile off
    profile on
end

[~, states, report] = simulateScheduleAD(initstate, model, schedule,...
                                                'OutputMinisteps', true,...
                                                'NonLinearSolver', nls);
if doprofiling
    profile viewer
end

%%  Process output

ind = cellfun(@(x) not(isempty(x)), states);
states = states(ind);
E = cellfun(@(x) x.(ctrl).E, states);
I = cellfun(@(x) x.(ctrl).I, states);
time = cellfun(@(x) x.time, states);

figure;
plot(time/hour, E);
xlabel('time  / h');
ylabel('voltage  / V');
title(modelcase);
grid on
axis tight
axis([0 max(time/hour) 3.75 4.15])

%% Plot an animated summary of the results
% plotDashboard(model, states, 'step', 0);

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
