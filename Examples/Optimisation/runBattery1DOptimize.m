%% Pseudo-Two-Dimensional (P2D) Lithium-Ion Battery Model
% This example demonstrates how to setup a P2D model of a Li-ion battery
% and run a simple simulation.

% clear the workspace and close open figures
clear
close all

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa optimization

%% Setup the properties of Li-ion battery materials and cell design

jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
jsonstruct.include_current_collectors = false;
jsonstruct.use_thermal = false;

% We define some shorthand names for simplicity.
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
am      = 'ActiveMaterial';
cc      = 'CurrentCollector';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
sep     = 'Separator';

jsonstruct.(ne).(am).diffusionModelType = 'simple';
jsonstruct.(pe).(am).diffusionModelType = 'simple';

inputparams = BatteryInputParams(jsonstruct);

inputparams.(ctrl).useCVswitch = true;

%% Setup the geometry and computational grid

gen = BatteryGeneratorP2D();

% Now, we update the inputparams with the properties of the grid.
inputparams = gen.updateBatteryInputParams(inputparams);

%%  Initialize the battery model.

model = Battery(inputparams);

model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', true);


%% Setup the time step schedule
% Smaller time steps are used to ramp up the current from zero to its
% operational value. Larger time steps are then used for the normal
% operation.

CRate = model.Control.CRate;

switch model.(ctrl).controlPolicy
  case 'CCCV'
    total = 3.5*hour/CRate;
  case 'CCDischarge'
    total = 1.2*hour/CRate;
  otherwise
    error('control policy not recognized');
end

n    = 40;
dt   = total*0.7/n;
step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

% Setup the control by assigning a source and stop function.

control = model.Control.setupScheduleControl();

nc = 1;
nst = numel(step.control);
ind = floor(((0 : nst - 1)/nst)*nc) + 1;

step.control = ind;
control.Imax = model.Control.Imax;
control = repmat(control, nc, 1);

schedule = struct('control', control, 'step', step);

%% Setup the initial state of the model
% The initial state of the model is dispatched using the
% model.setupInitialState() method.

initstate = model.setupInitialState();

%% Setup the properties of the nonlinear solver
nls = NonLinearSolver();
% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10;
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false;
% Change default tolerance for nonlinear solver
nls.timeStepSelector = SimpleTimeStepSelector();
model.nonlinearTolerance = 1e-3*model.Control.Imax;
% Set verbosity
model.verbose = true;

%% Run the simulation

[~, states, ~] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

%% Process output and recover the output voltage and current from the output states.

ind = cellfun(@(x) not(isempty(x)), states);
states = states(ind);

E    = cellfun(@(x) x.Control.E, states);
time = cellfun(@(x) x.time, states);

figure;
plot(time/hour, E, '*-', 'displayname', 'initial');
xlabel('time  / h');
ylabel('voltage  / V');
grid on

%% Calculate the energy
obj = @(model, states, schedule, varargin) EnergyOutput(model, states, schedule, varargin{:});
vals = obj(model, states, schedule);
totval = sum([vals{:}]);

% Compare with trapezoidal integral: they should be about the same
I = cellfun(@(x) x.Control.I, states);
totval_trapz = trapz(time, E.*I);
fprintf('Rectangle rule: %g Wh, trapezoidal rule: %g Wh\n', totval/hour, totval_trapz/hour);

%% Setup the optimization problem
state0 = initstate;
SimulatorSetup = struct('model', model, 'schedule', schedule, 'state0', state0);

parameters = {};

paramsetter = PorositySetter(model, {ne, sep, pe});

getporo = @(model, notused) paramsetter.getValues(model);
setporo = @(model, notused, v) paramsetter.setValues(model, v);

parameters = addParameter(parameters, SimulatorSetup, ...
                          'name'     , 'porosity', ...
                          'belongsTo', 'model'       , ...
                          'boxLims'  , [0.1, 0.9]    , ...
                          'location' , {''}          , ...
                          'getfun'   , getporo       , ...
                          'setfun'   , setporo);

setfun = @(x, location, v) struct('Imax', v                                                        , ...
                                  'src', @(time, I, E) rampupSwitchControl(time, 0.1, I, E, v, 3.0), ...
                                  'stopFunction', @(model,state,state_prev)(false)                 , ...
                                  'CCDischarge', true);

for i = 1 : max(schedule.step.control)
    parameters = addParameter(parameters, SimulatorSetup, ...
                              'name'        , 'Imax'                       , ...
                              'belongsTo'   , 'schedule'                   , ...
                              'boxLims'     , model.Control.Imax*[0.5, 2], ...
                              'location'    , {'control', 'Imax'}          , ...
                              'getfun'      , []                           , ...
                              'setfun'      , setfun                       , ...
                              'controlSteps', i);
end

objmatch = @(model, states, schedule, varargin) EnergyOutput(model, states, schedule, varargin{:});
fn       = @plotAfterStepIV;
obj      = @(p) evalObjectiveBattmo(p, objmatch, SimulatorSetup, parameters, 'objScaling', totval, 'afterStepFn', fn);

%%
doOptimization = true;
if doOptimization

    % The parameters must be scaled to [0,1]
    p_base = getScaledParameterVector(SimulatorSetup, parameters);
    p_base = p_base - 0.1;

    % Solve the optimization problem using BFGS. One can adjust the
    % tolerances and the maxIt option to see how it effects the
    % optimum.
    [v, p_opt, history] = unitBoxBFGS(p_base, obj, 'gradTol', 1e-7, 'objChangeTol', 1e-4, 'maxIt', 20);

    % Compute objective at optimum
    setup_opt = updateSetupFromScaledParameters(SimulatorSetup, parameters, p_opt);
    [~, states_opt, ~] = simulateScheduleAD(setup_opt.state0, setup_opt.model, setup_opt.schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);
    time_opt = cellfun(@(x) x.time, states_opt);
    E_opt = cellfun(@(x) x.Control.E, states_opt);
    I_opt = cellfun(@(x) x.Control.I, states_opt);
    totval_trapz_opt = trapz(time_opt, E_opt.*I_opt);

    % Print optimal parameters
    fprintf('Base and optimized parameters:\n');
    for k = 1:numel(parameters)
        % Get the original and optimized values
        p0 = parameters{k}.getParameter(SimulatorSetup);
        pu = parameters{k}.getParameter(setup_opt);

        % Print
        fprintf('%s\n', parameters{k}.name);
        fprintf('%g %g\n', p0, pu);
    end

    fprintf('Energy changed from %g to %g Wh\n', totval_trapz/hour, totval_trapz_opt/hour);

    % Plot
    figure; hold on; grid on
    E    = cellfun(@(x) x.Control.E, states);
    time = cellfun(@(x) x.time, states);
    plot(time/hour, E, '*-', 'displayname', 'initial');
    plot(time_opt/hour, E_opt, 'r*-', 'displayname', 'optimized');
    xlabel('time  / h');
    ylabel('voltage  / V');
    legend;
end


%%
doCompareGradient = true;
if doCompareGradient

    p = getScaledParameterVector(SimulatorSetup, parameters);
    [vad, gad]   = evalObjectiveBattmo(p, objmatch, SimulatorSetup, parameters, 'gradientMethod', 'AdjointAD');
    [vnum, gnum] = evalObjectiveBattmo(p, objmatch, SimulatorSetup, parameters, 'gradientMethod', 'PerturbationADNUM', 'PerturbationSize', 1e-5);

    fprintf('Gradient computed using adjoint:\n');
    display(gad);
    fprintf('Numerical gradient:\n');
    display(gnum);
    fprintf('Relative error:\n')
    display(abs(gnum-gad)./abs(gad));
end


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
