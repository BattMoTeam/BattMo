%% Pseudo-Two-Dimensional (P2D) Lithium-Ion Battery Model
% This example demonstrates how to setup a P2D model of a Li-ion battery
% and run a simple simulation.

% Clear the workspace and close open figures
clear
close all

%% Import the required modules from MRST
% load MRST modules

mrstModule add ad-core mrst-gui mpfa agmg linearsolvers

%% Setup the properties of Li-ion battery materials and cell design
% The properties and parameters of the battery cell, including the
% architecture and materials, are set using an instance of
% :class:`BatteryInputParams <Battery.BatteryInputParams>`. This class is
% used to initialize the simulation and it propagates all the parameters
% throughout the submodels. The input parameters can be set manually or
% provided in json format. All the parameters for the model are stored in
% the inputparams object.

jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));

% We define some shorthand names for simplicity.
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
co      = 'Coating';
am      = 'ActiveMaterial';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
cc      = 'CurrentCollector';

jsonstruct.use_thermal = false;
jsonstruct.include_current_collectors = false;

inputparams = BatteryInputParams(jsonstruct);

inputparams.(ne).(co).(am).diffusionModelType = 'simple';
inputparams.(pe).(co).(am).diffusionModelType = 'simple';

use_cccv = true;
if use_cccv
    cccvstruct = struct( 'controlPolicy'     , 'CCCV',  ...
                         'initialControl'    , 'discharging', ...
                         'CRate'             , 1.5          , ...
                         'DRate'             , 1            , ...
                         'lowerCutoffVoltage', 2.4          , ...
                         'upperCutoffVoltage', 4.1          , ...
                         'dIdtLimit'         , 0.01         , ...
                         'dEdtLimit'         , 0.01);
    cccvinputparams = CcCvControlModelInputParams(cccvstruct);
    inputparams.Control = cccvinputparams;
end


%% Setup the geometry and computational grid
% Here, we setup the 1D computational grid that will be used for the
% simulation. The required discretization parameters are already included
% in the class BatteryGeneratorP2D.
gen = BatteryGeneratorP2D();

% Now, we update the inputparams with the properties of the grid.
inputparams = gen.updateBatteryInputParams(inputparams);


%%  Initialize the battery model.
% The battery model is initialized by sending inputparams to the Battery class
% constructor. see :class:`Battery <Battery.Battery>`.
model = Battery(inputparams);

inspectgraph = false;
if inspectgraph
    cgti = model.computationalGraph;
    return
end

%% Compute the nominal cell capacity and choose a C-Rate
% The nominal capacity of the cell is calculated from the active materials.
% This value is then combined with the user-defined C-Rate to set the cell
% operational current.

DRate = model.Control.DRate;

%% Setup the time step schedule
% Smaller time steps are used to ramp up the current from zero to its
% operational value. Larger time steps are then used for the normal
% operation.
switch model.(ctrl).controlPolicy
  case 'CCCV'
    CRate = model.Control.CRate;
    total = 1.5*hour/(CRate + DRate);
  case 'CCDischarge'
    total = 1.4*hour/DRate;
  otherwise
    error('control policy not recognized');
end

n  = 50;
dt = total/n;
step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

% we setup the control by assigning a source and stop function.

control = model.Control.setupScheduleControl();

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step);

%% Setup the initial state of the model
% The initial state of the model is setup using the model.setupInitialState() method.

initstate = model.setupInitialState();

%% Setup the properties of the nonlinear solver
nls = NonLinearSolver();

linearsolver = 'amgcl';
switch linearsolver
  case 'amgcl'
    nls.LinearSolver = AMGCLSolverAD('verbose', true, 'reduceToCell', false);
    nls.LinearSolver.tolerance = 1e-4;
    nls.LinearSolver.maxIterations = 30;
    nls.maxIterations = 10;
    nls.verbose = 10;
  case 'battery'
    nls.LinearSolver = LinearSolverBatteryExtra('verbose'     , false, ...
                                                'reduceToCell', true, ...
                                                'verbosity'   , 3    , ...
                                                'reuse_setup' , false, ...
                                                'method'      , 'direct');
    nls.LinearSolver.tolerance = 1e-4;
  case 'direct'
    disp('standard direct solver')
  otherwise
    error('Unknown solver %s', linearsolver);
end

% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10;
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false;
nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);
% Change default tolerance for nonlinear solver
model.nonlinearTolerance = 1e-3*model.Control.Imax;
% Set verbosity
model.verbose = true;

%% Run the simulation
[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

%% Process output and recover the output voltage and current from the output states.
ind = cellfun(@(x) not(isempty(x)), states);
states = states(ind);
E = cellfun(@(x) x.Control.E, states);
I = cellfun(@(x) x.Control.I, states);
T = cellfun(@(x) max(x.(thermal).T), states);
Tmax = cellfun(@(x) max(x.ThermalModel.T), states);
% [SOCN, SOCP] =  cellfun(@(x) model.calculateSOC(x), states);
time = cellfun(@(x) x.time, states);

figure
plot(time/hour, E);
grid on
xlabel 'time  / h';
ylabel 'potential  / V';

%%
its  = getReportOutput(report, 'type', 'linearIterations');
nits = getReportOutput(report, 'type', 'nonlinearIterations');

figure(33),clf
plot(its.time, its.total./nits.total)
figure(44),clf,hold on
shift = 0;

for i = 1:numel(report.ControlstepReports)
    creport = report.ControlstepReports{i};
    for j = 1:numel(creport.StepReports)
        sreport = creport.StepReports{j};
        vits = [];
        for k = 1:numel(sreport.NonlinearReport)
            nreport = sreport.NonlinearReport{k};
            try
                its = nreport.LinearSolver.Iterations;
                its = nreport.LinearSolver.precondIterations_phi;
                vits = [vits, its];
            catch
            end
        end
        plot(shift+(1:numel(vits)), vits, '*-');
        shift = shift + numel(vits);
    end
end

figure(55)
plot(nits.total);

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
