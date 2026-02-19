clear all
close all


%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa agmg linearsolvers


jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));

% We define some shorthand names for simplicity.
elyte   = 'Electrolyte';
am      = 'ActiveMaterial';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';

jsonstruct.use_thermal = false;
jsonstruct.include_current_collectors = false;

jsonstruct.(pe).(am).diffusionModelType = 'full';
jsonstruct.(ne).(am).diffusionModelType = 'full';

inputparams = BatteryInputParams(jsonstruct);

gen = BatteryGenerator1D();

% Now, we update the inputparams with the properties of the mesh. 
inputparams = gen.updateBatteryInputParams(inputparams);
G = inputparams.(ne).(am).G;

batterymodel = Battery(inputparams);
batterymodel = batterymodel.setupComputationalGraph();
cgtbattery   = batterymodel.computationalGraph;

hcjsonstruct.(am) = jsonstruct.(ne).(am);
hcjsonstruct.(elyte) = jsonstruct.Electrolyte;
hcjsonstruct.(ctrl) = jsonstruct.(ctrl);


inputparams = HalfCellInputParams(hcjsonstruct);
inputparams.(am).G    = G;

model = HalfCell(inputparams);

%% Defining the SOC
model.SOC = 0.99;

modelHC = model.setupComputationalGraph();
cgt = modelHC.computationalGraph;
[g, edgelabels] = getComputationalGraph(cgt);
plot(g)

modelHC.AutoDiffBackend= AutoDiffBackend();

inspectgraph = false;
if inspectgraph
    % plot the computational graph
    cgt = ComputationalGraphTool(modelHC);
    cgt.getComputationalGraph('doplot', true);
    return
end

%% Compute the nominal cell capacity and choose a C-Rate
% The nominal capacity of the cell is calculated from the active materials.
% This value is then combined with the user-defined C-Rate to set the cell
% operational current. 

CRate = model.Control.CRate;

%% Setup the time step schedule 
% Smaller time steps are used to ramp up the current from zero to its
% operational value. Larger time steps are then used for the normal
% operation.
switch model.(ctrl).controlPolicy
  case 'CCCV'
    total = 3.5*hour/CRate;
  case 'IEswitch'
    total = 1.4*hour/CRate;
  otherwise
    error('control policy not recognized');
end

n  = 100;
dt = total/n;
step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

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
schedule = struct('control', control, 'step', step); 

%% Setup the initial state of the model
% The initial state of the model is setup using the model.setupInitialState() method.

initstate = model.setupInitialState(); 

%% Setup the properties of the nonlinear solver 
nls = NonLinearSolver();

linearsolver = 'direct';
switch linearsolver
  case 'agmg'
    mrstModule add agmg
    nls.LinearSolver = AGMGSolverAD('verbose', true, 'reduceToCell', false); 
    nls.LinearSolver.tolerance = 1e-3; 
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
    error()
end

% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10;
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false;
nls.timeStepSelector=StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);
% Change default tolerance for nonlinear solver
model.nonlinearTolerance = 1e-3*model.Control.Imax;
% Set verbosity
model.verbose = true;

%% Run the simulation
[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 




%% Plotting
ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
E = cellfun(@(x) x.Control.E, states); 
I = cellfun(@(x) x.Control.I, states);
time = cellfun(@(x) x.time, states);

plot(time,E);



