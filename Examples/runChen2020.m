%% Battery 1D model
% Include presentation of the test case (use rst format)

% load MRST modules
mrstModule add ad-core multimodel mrst-gui battery mpfa

% We create an instance of BatteryInputParams. This class is used to initiate the battery simulator and it propagates
% all the parameters through out the submodels.

% The input parameters can be given in json format. The json file is read and used to populate the paramobj object.
jsonstruct = parseBatmoJson('JsonDatas/Chen2020/chenBattery.json');

paramobj = BareBatteryInputParams(jsonstruct);

% Some shortcuts used for the sub-models
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
elyte   = 'Electrolyte';

%% We setup the battery geometry.
% Here, we use a 1D model and the class BatteryGenerator1D already contains the discretization parameters
gen = BareBatteryGenerator1D();
% We update pamobj with grid data
paramobj = gen.updateBatteryInputParams(paramobj);

%%  The Battery model is initialized by sending paramobj to the Battery class constructor 

model = BareBattery(paramobj);

%% We compute the cell capacity and chose a discharge rate
C      = computeCellCapacity(model);
CRate  = 1/5; 
inputI = (C/hour)*CRate; % current 

%% We setup the schedule 
% We use different time step for the activation phase (small time steps) and the following discharging phase

% We start with rampup time steps to go through the activation phase 
dt1   = rampupTimesteps(0.1, 0.1, 10);
% We choose time steps for the rest of the simulation (discharge phase)
dt2   = 0.1*hour*ones(30, 1);
% We concatenate the time steps
dt    = [dt1; dt2];
times = [0; cumsum(dt)]; 
tt    = times(2 : end); 
step  = struct('val', diff(times), 'control', ones(numel(tt), 1)); 

% stopFunc = @(model, state, state_prev) (state.(pe).E < 2.0); 
stopFunc = @(model, state, state_prev) (false);

tup = 0.1; % rampup value for the current function, see rampupSwitchControl
inputE = 3.6; % Value when current control switches to voltage control
srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, inputI, inputE);

% we setup the control by assigning a source and stop function.
control = repmat(struct('src', srcfunc, 'stopFunction', stopFunc), 1, 1); 

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step); 

%%  We setup the initial state
initstate = model.setupInitialState(); 

% Setup nonlinear solver 
nls = NonLinearSolver(); 
% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10; 
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false; 
% Change default tolerance for nonlinear solver
model.nonlinearTolerance = 1e-5; 
% Set verbosity
model.verbose = false;

% Run simulation

[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 