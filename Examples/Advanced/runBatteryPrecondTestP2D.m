%% Pseudo-Two-Dimensional (P2D) Lithium-Ion Battery Model
% This example demonstrates how to setup a P2D model of a Li-ion battery
% and run a simple simulation.

if mrstPlatform('octave')
    error('This demo cannot be run from Octave since it does not support the linear solvers');
end

% Clear the workspace and close open figures
clear
close all

mrstDebug(20);

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
am      = 'ActiveMaterial';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
cc      = 'CurrentCollector';

jsonstruct.use_thermal = true;
jsonstruct.include_current_collectors = false;

diffusionModelType = 'full';

jsonstruct.(pe).(am).diffusionModelType = diffusionModelType;
jsonstruct.(ne).(am).diffusionModelType = diffusionModelType;

inputparams = BatteryInputParams(jsonstruct);

% inputparams.(ne).(am).(sd).N = 5;
% inputparams.(pe).(am).(sd).N = 5;

use_cccv = false;
if use_cccv
    cccvstruct = struct( 'controlPolicy'     , 'CCCV',  ...
                         'CRate'             , 1.5 , ...
                         'DRate'             , 1   , ...
                         'lowerCutoffVoltage', 2.4 , ...
                         'upperCutoffVoltage', 4.1 , ...
                         'dIdtLimit'         , 0.01, ...
                         'dEdtLimit'         , 0.01);
    cccvinputparams = CcCvControlModelInputParams(cccvstruct);
    inputparams.Control = cccvinputparams;
end


%% Setup the geometry and computational grid
% Here, we setup the 1D computational grid that will be used for the
% simulation. The required discretization parameters are already included
% in the class BatteryGeneratorP2D.
gen = BatteryGeneratorP2D();
gen.resolutionFactor = 100;
gen = gen.applyResolutionFactors();

% Now, we update the inputparams with the properties of the grid.
inputparams = gen.updateBatteryInputParams(inputparams);


%%  Initialize the battery model.
% The battery model is initialized by sending inputparams to the Battery class
% constructor. see :class:`Battery <Battery.Battery>`.
model = Battery(inputparams);

inspectgraph = false;
if inspectgraph
    % plot the computational graph
    cgti = ComputationalGraphInteractiveTool(model);
    cgti.getComputationalGraph('doplot', true);
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
    total = 1.5*hour/(DRate + CRate);
  case 'CCDischarge'
    total = 1.4*hour/DRate;
  otherwise
    error('control policy not recognized');
end

n  = 100;
dt = total/n;

step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

% Setup the control by assigning a source and stop function.
control = model.Control.setupScheduleControl();

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step);

%% Setup the initial state of the model
% The initial state of the model is setup using the model.setupInitialState() method.

initstate = model.setupInitialState();

%% Setup the properties of the nonlinear solver
nls = NonLinearSolver();

clear setup

casenumber = 4;

switch casenumber

  case 1

    setup.library = 'matlab';
    setup.method = 'direct';

  case 2

    setup.method = 'gmres';
    setup.options.method = 'grouped';
    setup.options.solverspec.name = 'agmg'; % not used for now

  case 3

    setup.method = 'gmres';
    setup.options.method = 'separate';
    solvers = {};
    solverspec.variable = 'phi';
    solverspec.solverspec.name = 'direct';
    solvers{end + 1} = solverspec;
    solverspec.variable = 'c';
    solverspec.solverspec.name = 'direct';
    solvers{end + 1} = solverspec;
    solverspec.variable = 'T';
    solverspec.solverspec.name = 'direct';
    solvers{end + 1} = solverspec;
    setup.options.solvers = solvers;

  case 4

    switch diffusionModelType
      case 'simple'
        jsonfilename = fullfile(battmoDir, 'Utilities/JsonSchemas/Tests/linearsolver3.json');
      case 'full'
        jsonfilename = fullfile(battmoDir, 'Utilities/JsonSchemas/Tests/linearsolver4.json');
      otherwise
        error('diffusionModelType not covered')
    end
    jsonsrc = fileread(jsonfilename);
    setup = battMojsondecode(jsonsrc);

  otherwise

    error('case number not recognized');

end

if isfield(setup, 'reduction')
    model = model.setupSelectedModel('reduction', setup.reduction);
end

nls.LinearSolver = BatteryLinearSolver('verbose'          , 2    , ...
                                       'reuse_setup'      , false, ...
                                       'linearSolverSetup', setup);

if isfield(nls.LinearSolver.linearSolverSetup, 'gmres_options')
    nls.LinearSolver.linearSolverSetup.gmres_options.tol = 1e-3*model.Control.Imax;
    nls.LinearSolver.linearSolverSetup.gmres_options.maxit = 10;
end

% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 20;
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false;
% nls.timeStepSelector=StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);
% Change default tolerance for nonlinear solver
model.nonlinearTolerance = 1e-3*model.Control.Imax;
% Set verbosity
model.verbose = true;

%% Run the simulation
% profile off
% profile on
[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);
% profile off
% profile viewer

%% Process output and recover the output voltage and current from the output states.
ind = cellfun(@(x) not(isempty(x)), states);
states = states(ind);
E = cellfun(@(x) x.Control.E, states);
I = cellfun(@(x) x.Control.I, states);
Tmax = cellfun(@(x) max(x.ThermalModel.T), states);
% [SOCN, SOCP] =  cellfun(@(x) model.calculateSOC(x), states);
time = cellfun(@(x) x.time, states);

figure
plot(time, E);


%%

its = getReportOutput(report,'type','linearIterations')
nits= getReportOutput(report,'type','nonlinearIterations')

figure
plot(its.time, its.total./nits.total)
xlabel('time');
title('linear iterations/newton iterations')


switch setup.method

  case "grouped-gmres"

    pits = zeros(numel(its.time), 1);

    counter = 1;

    for icontrol = 1 : numel(report.ControlstepReports)

        creport = report.ControlstepReports{icontrol};

        for istep = 1 : numel(creport.StepReports)

            sreport = creport.StepReports{istep};

            for k = 1 : numel(sreport.NonlinearReport);

                nreport = sreport.NonlinearReport{k};
                try
                    it = nreport.LinearSolver.precondReports.Iterations;
                    pits(counter) = pits(counter) + it;
                end

            end
        end
    end

    figure
    plot(its.time, pits./its.total);
    xlabel('time');
    title('preconditioner iteration/linear iteration');

  case "separate-variable-gmres"

    npreconds = numel(setup.preconditioners);
    pits = zeros(numel(its.time), npreconds);

    counter = 1;

    for icontrol = 1 : numel(report.ControlstepReports)

        creport = report.ControlstepReports{icontrol};

        for istep = 1 : numel(creport.StepReports)

            sreport = creport.StepReports{istep};

            for k = 1:numel(sreport.NonlinearReport);

                nreport = sreport.NonlinearReport{k};

                for iprecond = 1 : npreconds

                    try
                        it = nreport.LinearSolver.precondReports{iprecond}.Iterations;
                        pits(counter, iprecond) = pits(counter, iprecond) + it;
                    end

                end
            end

            counter = counter + 1;
        end
    end

    for iprecond = 1 : npreconds

        figure
        plot(its.time, pits(:, iprecond)./its.total);
        xlabel('time');
        switch iprecond
          case 1
            titlestr = 'phi';
          case 2
            titlestr = 'c';
          case 3
            titlestr = 'T';
        end
        titlestr = sprintf('preconditioner iteration/linear iteration - %s', titlestr);
        title(titlestr);

    end

end



%% Plot the the output voltage and current
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
