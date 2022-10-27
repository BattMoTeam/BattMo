%% Pseudo-Two-Dimensional (P2D) Lithium-Ion Battery Model
% This example demonstrates how to setup a P2D model of a Li-ion battery
% and run a simple simulation.

% clear the workspace and close open figures
clear
close all
clc


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

paramobj = BatteryInputParams(jsonstruct);

%% Setup the geometry and computational mesh

gen = BatteryGenerator1D();

% Now, we update the paramobj with the properties of the mesh. 
paramobj = gen.updateBatteryInputParams(paramobj);

paramobj.(ne).(cc).EffectiveElectricalConductivity = 1e5;
paramobj.(pe).(cc).EffectiveElectricalConductivity = 1e5;

%%  Initialize the battery model. 

model = Battery(paramobj);
model.AutoDiffBackend= AutoDiffBackend();

%% Compute the nominal cell capacity and choose a C-Rate

CRate = model.Control.CRate;

%% Setup the time step schedule 
% Smaller time steps are used to ramp up the current from zero to its
% operational value. Larger time steps are then used for the normal
% operation.
switch model.(ctrl).controlPolicy
  case 'CCCV'
    total = 3.5*hour/CRate;
  case 'IEswitch'
    total = 1.2*hour/CRate;
  otherwise
    error('control policy not recognized');
end

n     = 40;
dt    = total*0.7/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

% we setup the control by assigning a source and stop function.
% control = struct('CCCV', true); 
%  !!! Change this to an entry in the JSON with better variable names !!!

switch model.Control.controlPolicy
  case 'IEswitch'
    tup = 0.1; % rampup value for the current function, see rampupSwitchControl
    srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                model.Control.Imax, ...
                                                model.Control.lowerCutoffVoltage);
    % we setup the control by assigning a source and stop function.
    control = struct('src', srcfunc, 'IEswitch', true);
  case 'CCCV'
    control = struct('CCCV', true);
  otherwise
    error('control policy not recognized');
end

% This control is used to set up the schedule

nc = 1;
nst = numel(step.control);
ind = floor(([0 : nst - 1]/nst)*nc) + 1;

step.control = ind;
control.Imax = model.Control.Imax;
control = repmat(control,nc,1);
schedule = struct('control', control, 'step', step); 

%% Setup the initial state of the model
% The initial state of the model is dispatched using the
% model.setupInitialState()method. 

initstate = model.setupInitialState(); 

%% Setup the properties of the nonlinear solver 
nls = NonLinearSolver(); 
% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10; 
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false; 
%nls.timeStepSelector=StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);
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

plot(time, E, '*-');

%% Plot the the output voltage and current
% gradients to controlmodel.use_thermal =false;

obj = @(model, states, schedule, varargin) EnergyOutput(model, states, schedule,varargin{:});%,'step',step);
vals = obj(model, states, schedule);
totval = sum([vals{:}]);

% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 30 iterations

state0 = initstate;
SimulatorSetup = struct('model', model, 'schedule', schedule, 'state0', state0);

parameters={};

part = ones(model.(elyte).(sep).G.cells.num, 1);
parameters{end+1} = ModelParameter(SimulatorSetup, ...
                                   'name'     , 'porosity_sep'           , ...
                                   'belongsTo', 'model'                  , ...
                                   'boxLims'  , [0.1 0.9]                , ...
                                   'lumping'  , part                     , ...
                                   'location' , {elyte, sep, 'porosity'} , ...
                                   'getfun'   , []                       , ...
                                   'setfun'   , []);

part = ones(model.(ne).(am).G.cells.num, 1);
parameters{end+1} = ModelParameter(SimulatorSetup, ...
                                   'name'     , 'porosity_nam'      , ...
                                   'belongsTo', 'model'             , ...
                                   'boxLims'  , [0.07 0.5]          , ...
                                   'lumping'  , part                , ...
                                   'location' , {ne, am, 'porosity'}, ...
                                   'getfun'   , []                  , ...
                                   'setfun'   , []);

part = ones(model.(pe).(am).G.cells.num, 1);
parameters{end+1} = ModelParameter(SimulatorSetup, ...
                                   'name'     , 'porosity_pam'      , ...
                                   'belongsTo', 'model'             , ...
                                   'boxLims'  , [0.07 0.5]          , ...
                                   'lumping'  , part                , ...
                                   'location' , {pe, am,'porosity'} , ...
                                   'getfun'   , []                  , ...
                                   'setfun'   , []);

setfun = @(x, location, v) struct('Imax', v, ...
                                  'src', @(time, I, E) rampupSwitchControl(time, 0.1, I, E, v, 3.0), ...
                                  'IEswitch', true);

boxlims = model.Control.Imax*[0.5, 1.5];

for i = 1 : max(schedule.step.control)
    parameters{end+1} = ModelParameter(SimulatorSetup, ...
                                       'name'        , 'Imax'             , ...
                                       'belongsTo'   , 'schedule'         , ...
                                       'boxLims'     , boxlims            , ...
                                       'location'    , {'control', 'Imax'}, ...
                                       'getfun'      , []                 , ...
                                       'setfun'      , setfun             , ...
                                       'controlSteps', i);
end

objmatch = @(model, states, schedule, varargin) EnergyOutput(model, states, schedule, varargin{:});
fn       = @plotAfterStepIV;
obj      = @(p) evalMatchBattmo(p, objmatch, SimulatorSetup, parameters, 'objScaling', -totval, 'afterStepFn', fn);

doOptimization = true;
if doOptimization
    
    p_base = getScaledParameterVector(SimulatorSetup, parameters);
    p_base = p_base - 0.1;
    [v, p_opt, history] = unitBoxBFGS(p_base, obj, 'gradTol', 1e-7, 'objChangeTol', 1e-4);
    
end


doCompareGradient = false;
if doCompareGradient
    
    p = getScaledParameterVector(SimulatorSetup, parameters);
    [vad, gad]   = evalMatchBattmo(p, objmatch, SimulatorSetup, parameters, 'Gradient', 'AdjointAD');
    [vnum, gnum] = evalMatchBattmo(p, objmatch, SimulatorSetup, parameters, 'Gradient', 'PerturbationADNUM', 'PerturbationSize', 1e-5);

    fprintf('Gradient computed using adjoint:\n');
    display(gad);
    fprintf('Numerical gradient:\n');
    display(gnum);
    
end