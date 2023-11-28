if mrstPlatform('octave')
    error('This demo cannot be run from Octave since it does not support the linear solvers');
end

clear all
close all

% mrstDebug(20);

% Setup mrst modules
mrstModule add ad-core mrst-gui mpfa agmg linearsolvers

%% We setup the geometrical parameters for a 4680 battery.
%% Those will be gathered in structure spiralparams (see below) and used by SpiralBatteryGenerator to generate the spiral layered geometry of the jelly roll

% Inner radius of the jelly roll
rInner = 2*milli*meter;

% widths of each component ordered as
% - positive current collector
% - positive electrode
% - electrolyte separator
% - negative electrode
% - negative current collector

widths = [25, 64, 15, 57, 15]*micro*meter;

widthDict = containers.Map(...
    {'Separator'               , ...
     'NegativeCoating'         , ...
     'NegativeCurrentCollector', ...
     'PositiveCoating'         , ...
     'PositiveCurrentCollector'}, ...
    widths);

nwidths = [widthDict('PositiveCoating')          ;...
           widthDict('PositiveCurrentCollector') ;...
           widthDict('PositiveCoating')          ;...
           widthDict('Separator')                ;...
           widthDict('NegativeCoating')          ;...
           widthDict('NegativeCurrentCollector') ;...
           widthDict('NegativeCoating')          ;...
           widthDict('Separator')];

dr = sum(nwidths);

% Outer radius of the jelly roll
rOuter = 22.36*milli*meter;
% Height of the jelly roll
L = 80*milli*meter;

dR = rOuter - rInner;
% Computed number of windings
nwindings = ceil(dR/dr);

% number of discretization cells in radial direction for each component.
nrDict = containers.Map({'Separator'               , ...
                         'NegativeCoating'         , ...
                         'NegativeCurrentCollector', ...
                         'PositiveCoating'         , ...
                         'PositiveCurrentCollector'}, ...
                        [3, 3, 3, 3, 3]);

% Number of discretization cells in the angular direction
nas = 10;

% Number of discretization cells in the longitudonal
nL = 3;

% structure that describes the tab setups (see SpiralBatteryGenerator)
tabparams.tabcase   = 'aligned tabs';
tabparams.width     = 3*milli*meter;
tabparams.fractions = linspace(0.01, 0.9, 6);

testing = true;
if testing
    fprintf('We setup a smaller case for quicker testing\n');
    rOuter = rInner + 1*milli*meter;
dR = rOuter - rInner;
% Computed number of windings
nwindings = ceil(dR/dr);
nL = 2;
end

spiralparams = struct('nwindings'   , nwindings, ...
                      'rInner'      , rInner   , ...
                      'widthDict'   , widthDict, ...
                      'nrDict'      , nrDict   , ...
                      'nas'         , nas      , ...
                      'L'           , L        , ...
                      'nL'          , nL       , ...
                      'tabparams'   , tabparams, ...
                      'angleuniform', true);

% The input material parameters given in json format are used to populate the paramobj object.
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

jsonstruct.include_current_collectors = true;

jsonstruct.use_thermal = true;

diffusionModelType = 'full';

jsonstruct.(pe).(co).(am).diffusionModelType = diffusionModelType;
jsonstruct.(ne).(co).(am).diffusionModelType = diffusionModelType;

paramobj = BatteryInputParams(jsonstruct);

% paramobj.(ne).(am).(sd).N = 5;
% paramobj.(pe).(am).(sd).N = 5;

paramobj.(ctrl).lowerCutoffVoltage = 3;
paramobj.(ctrl).CRate = 0.1;

% th = 'ThermalModel';
% paramobj.(th).externalHeatTransferCoefficientSideFaces = 100*watt/meter^2;
% paramobj.(th).externalHeatTransferCoefficientTopFaces = 10*watt/meter^2;

gen = SpiralBatteryGenerator();

paramobj = gen.updateBatteryInputParams(paramobj, spiralparams);

model = Battery(paramobj);
model.AutoDiffBackend= AutoDiffBackend();

%% Setup schedule

CRate = model.(ctrl).CRate;
fac   = 2;
total = 1.1*hour/CRate;
n     = 10;
dt0   = total*1e-6;
times = getTimeSteps(dt0, n, total, fac);

% times = times(times < 4200);

%% We compute the cell capacity, which used to compute schedule from CRate

step = struct('val', diff(times), 'control', ones(numel(times) - 1, 1));

tup = 0.1/CRate;

simcase = 'discharge';

switch simcase

  case 'discharge'

    srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                model.Control.Imax, ...
                                                model.Control.lowerCutoffVoltage);
    % we setup the control by assigning a source and stop function.
    control = struct('src', srcfunc, 'CCDischarge', true);
    schedule  = struct('control', control, 'step', step);

    %% We setup the initial state
    initstate = model.setupInitialState();
    initElytec = 1*mol/litre;
    initstate.Electrolyte.c = initElytec*ones(model.Electrolyte.G.cells.num, 1);

  case 'charge'

    model.SOC = 0.01;
    initstate = model.setupInitialState();
    srcfunc  = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                 - model.Control.Imax, ...
                                                 model.Control.lowerCutoffVoltage);
    control = struct('src', srcfunc, 'CCDischarge', true);
    schedule = struct('control', control, 'step', step);

  otherwise
    error('simcase not recognized')

end

%% Setup the properties of the nonlinear solver
nls = NonLinearSolver();

clear setup

casenumber = 3;
% casenumber = 3;

switch casenumber

  case 1

    setup.library = 'matlab';

  case 2

    setup.method = 'gmres';
    setup.options.method = 'grouped';
    setup.options.solverspec.name = 'agmg'; % not used for now

  case 3

    switch diffusionModelType
      case 'simple'
        % jsonfilename = fullfile(battmoDir, 'Utilities/JsonSchemas/Tests/linearsolver3.json');
        jsonfilename = fullfile(battmoDir, 'Utilities/JsonSchemas/Tests/linearsolver5.json');
      case 'full'
        jsonfilename = fullfile(battmoDir, 'Utilities/JsonSchemas/Tests/linearsolver4.json');
      otherwise
        error('diffusionModelType not covered')
    end
    jsonsrc = fileread(jsonfilename);
    setup = jsondecode(jsonsrc);


  case 4

    switch diffusionModelType
      case 'simple'
        % jsonfilename = fullfile(battmoDir, 'Utilities/JsonSchemas/Tests/linearsolver3.json');
        jsonfilename = fullfile(battmoDir, 'Utilities/JsonSchemas/Tests/linearsolver5.json');
      case 'full'
        jsonfilename = fullfile(battmoDir, 'Utilities/JsonSchemas/Tests/linearsolver4.json');
      otherwise
        error('diffusionModelType not covered')
    end
    jsonsrc = fileread(jsonfilename);
    setup = jsondecode(jsonsrc);

    clear solver
    solver.library = 'matlab';
    solver.method  = 'direct';
    for iprec = 1 : numel(setup.preconditioners)
        setup.preconditioners(iprec).solver =  solver;
    end

    setup.verbose = 0;

  otherwise

    error('case number not recognized');

end

if isfield(setup, 'reduction')
    model = model.setupSelectedModel('reduction', setup.reduction);
end

if testing & (casenumber == 3)
    % We remove verbosity in case of small model
    setup.verbose = 0;
    if isfield(setup, 'preconditioners')
        pcs = setup.preconditioners;
        for ipc = 1 : numel(pcs)
            pc = pcs(ipc);
            pc.solver.verbose = 0;
            pc.solver.solver.verbose = false;
            pcs(ipc) = pc;
        end
        setup.preconditioners = pcs;
    end
end

%%
nls.LinearSolver = BatteryLinearSolver('verbose'          , 0    , ...
                                       'reuse_setup'      , false, ...
                                       'linearSolverSetup', setup);

if isfield(nls.LinearSolver.linearSolverSetup, 'gmres_options')
    nls.LinearSolver.linearSolverSetup.gmres_options.tol = 1e-3*model.Control.Imax;
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

dopacked = true;

if dopacked
    % Run simulation
    dataFolder = 'BattMo';
    problem = packSimulationProblem(initstate, model, schedule, dataFolder, ...
                                    'Name'           , 'jellyroll', ...
                                    'NonLinearSolver', nls);
    problem.SimulatorSetup.OutputMinisteps = true;

    clearSimulation = true;
    if clearSimulation
        %% clear previously computed simulation
        clearPackedSimulatorOutput(problem, 'prompt', false);
    end
    simulatePackedProblem(problem);
    [globvars, states, reports] = getPackedSimulatorOutput(problem);

else

    fn = afterStepConvergencePlots(nls);

    [~, states, reports] = simulateScheduleAD(initstate, model, schedule, ...
                                              'NonLinearSolver', nls, ...
                                              'afterStepFn'    , fn);

    reports = reports.ControlstepReports;

end

%% Process output and recover the output voltage and current from the output states.

ind = cellfun(@(x) not(isempty(x)), states);
states = states(ind);
E = cellfun(@(x) x.Control.E, states);
I = cellfun(@(x) x.Control.I, states);
time = cellfun(@(x) x.time, states);

figure
plot(time, E, 'linewidth', 3);
set(gca, 'fontsize', 18);
title('Cell Voltage / V')
xlabel('time')




%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
