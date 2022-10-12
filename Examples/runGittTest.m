%% Example: Galvanostatic Intermittent Charge Titration Technique (GITT) Simulation of a Lithium-Ion Battery Model
% This example demonstrates how to apply a GITT protocol in a Li-ion
% battery simulation.

% clear the workspace and close open figures
clear all
close all
clc


%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa
mrstVerbose off

%% Setup the properties of Li-ion battery materials and cell design
% The properties and parameters of the battery cell, including the
% architecture and materials, are set using an instance of
% :class:`BatteryInputParams <Battery.BatteryInputParams>`. This class is
% used to initialize the simulation and it propagates all the parameters
% throughout the submodels. The input parameters can be set manually or
% provided in json format. All the parameters for the model are stored in
% the paramobj object.
jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
paramobj = BatteryInputParams(jsonstruct);

% We define some shorthand names for simplicity.
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
eac     = 'ActiveMaterial';
cc      = 'CurrentCollector';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
ctrl    = 'Control';

%% Setup the geometry and computational mesh
% Here, we setup the geometry and computational mesh that will be used for
% the simulation. The user can select the dimensionality of the model. The
% required discretization parameters are already included in the associated
% battery generator class.
modelcase = '2D';

% Generate the battery based on the selected dimensionality in modelcase
switch modelcase
  case '1D'

    gen = BatteryGenerator1D();
    paramobj = gen.updateBatteryInputParams(paramobj);
    paramobj.(ne).(cc).EffectiveElectricalConductivity = 100;
    paramobj.(pe).(cc).EffectiveElectricalConductivity = 100;

    paramobj.(thermal).externalHeatTransferCoefficient = 1000;
    paramobj.(thermal).externalTemperature = paramobj.initT;

  case '2D'
    gen = BatteryGenerator2D();
    paramobj = gen.updateBatteryInputParams(paramobj);

    paramobj.(ne).(cc).EffectiveElectricalConductivity = 1e5;
    paramobj.(pe).(cc).EffectiveElectricalConductivity = 1e5;

    paramobj.(thermal).externalTemperature = paramobj.initT;
    paramobj.SOC = 0.99;

  case '3D'
    gen = BatteryGenerator3D();

    fac = 1;
    gen.facx = fac;
    gen.facy = fac;
    gen.facz = fac;
    gen = gen.applyResolutionFactors();
    paramobj = gen.updateBatteryInputParams(paramobj);

    paramobj.(thermal).externalTemperature = paramobj.initT;

end

%%  Initialize the battery model.
% The battery model is initialized by sending paramobj to the Battery class
% constructor. see :class:`Battery <Battery.Battery>`.
model = Battery(paramobj);

%% Compute the nominal cell capacity and choose a C-Rate
% The nominal capacity of the cell is calculated from the active materials.
% This value is then combined with the user-defined C-Rate to set the cell
% operational current.
C      = computeCellCapacity(model);
CRate  = 2;
inputI = (C/hour)*CRate;
inputE = model.Control.upperCutoffVoltage;

%% Setup the parameters of the GITT protocol
pulseFraction  = 0.01;
relaxationTime = 4*hour;
switchTime     = 1*milli*second; % switching time (linear interpolation between the two states)
dischargeTime  = pulseFraction*CRate*hour; % time of discharging

% Discretization parameters
numberOfIntervals             = 1/pulseFraction;
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
tup          = 1*milli*second; % rampup time
srcfunc_init = @(time, I, E) rampupSwitchControl(time, tup, I, E, inputI, inputE);
srcfunc_gitt = @(time, I, E) tabulatedIControl(time, tpoints, Ipoints);

control(1) = struct('src', srcfunc_init, 'IEswitch', true);
control(2) = struct('src', srcfunc_gitt, 'IEswitch', true);

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
use_diagonal_ad = false;
if(use_diagonal_ad)
    model.AutoDiffBackend = DiagonalAutoDiffBackend();
    model.AutoDiffBackend.useMex = true;
    model.AutoDiffBackend.modifyOperators = true;
    model.AutoDiffBackend.rowMajor = true;
    model.AutoDiffBackend.deferredAssembly = false; % error with true for now
end

use_iterative = false;
if(use_iterative)
    % nls.LinearSolver = LinearSolverBattery('method', 'iterative');
    % nls.LinearSolver = LinearSolverBattery('method', 'direct');
    mrstModule add agmg
    nls.LinearSolver = LinearSolverBattery('method', 'agmg', 'verbosity', 1);
    nls.LinearSolver.tol = 1e-3;
    nls.verbose = 10
end
model.nonlinearTolerance = 1e-5;
model.verbose = false;

%% Run simulation
doprofiling = false;
if doprofiling
    profile off
    profile on
end

[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule,...
                                                'OutputMinisteps', true,...
                                                'NonLinearSolver', nls);

%%  Process output

ind = cellfun(@(x) not(isempty(x)), states);
states = states(ind);
E = cellfun(@(x) x.(ctrl).E, states);
I = cellfun(@(x) x.(ctrl).I, states);
time = cellfun(@(x) x.time, states);

%% Plot an animated summary of the results
plotDashboard(model, states, 'step', 0);

%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
