
% clear the workspace and close open figures
clear all
close all

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa matlab_bgl

%% shortcuts

ne = 'NegativeElectrode';
pe = 'PositiveElectrode';
co = 'Coating';

am1 = 'ActiveMaterial1';
am2 = 'ActiveMaterial2';

bd = 'Binder';
ad = 'ConductingAdditive';

sd  = 'SolidDiffusion';
itf = 'Interface';

%% Setup the properties of Li-ion battery materials and cell design
jsonstruct_composite_material = parseBattmoJson('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_silicon_graphite.json');
jsonstruct_cell               = parseBattmoJson('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json');

jsonstruct_cell.(ne).(co) = rmfield(jsonstruct_cell.(ne).(co), 'ActiveMaterial');

jsonstruct = mergeJsonStructs({jsonstruct_composite_material, ...
                               jsonstruct_cell});

jsonstruct.use_thermal = false;
jsonstruct.include_current_collectors = false;

paramobj = BatteryInputParams(jsonstruct);

paramobj.(ne).(co).(am1).massFraction = 0.9;
paramobj.(ne).(co).(am2).massFraction = 0.08;
paramobj.(ne).(co).(bd).massFraction  = 0.01;
paramobj.(ne).(co).(ad).massFraction  = 0.01;

paramobj.Control.CRate = 0.1;

Paramobj = paramobj.validateInputParams();

gen = BatteryGenerator1D();

gen.xlength(4) = 1.8619*gen.xlength(4);

% gen.fac = 100;
% gen = gen.applyResolutionFactors();

% Now, we update the paramobj with the properties of the mesh. 
paramobj = gen.updateBatteryInputParams(paramobj);

model = Battery(paramobj);

model = model.validateModel();

model.AutoDiffBackend= AutoDiffBackend();

inspectgraph = false;
if inspectgraph
    cgt = model.computationalGraph;
    return
end

% model.Control.Imax = 1e-1;

%% Setup the initial state of the model
% The initial state of the model is setup using the model.setupInitialState() method.

initstate = model.setupInitialState();

%% setup schedule

%% Compute the nominal cell capacity and choose a C-Rate
% The nominal capacity of the cell is calculated from the active materials.
% This value is then combined with the user-defined C-Rate to set the cell
% operational current. 

CRate = model.Control.CRate;

%% Setup the time step schedule 
% Smaller time steps are used to ramp up the current from zero to its
% operational value. Larger time steps are then used for the normal
% operation.

ctrl = 'Control';

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

%% Run simulation

initstate = model.setupInitialState(); 

%% Setup the properties of the nonlinear solver 
nls = NonLinearSolver();

CRate = model.Control.CRate;

ctrl = 'Control';

switch model.(ctrl).controlPolicy
  case 'CCCV'
    total = 3.5*hour/CRate;
  case 'IEswitch'
    total = 1.4*hour/CRate;
  otherwise
    error('control policy not recognized');
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
tic
[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 
toc

%% plotting

set(0, 'defaultlinelinewidth', 3);

ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
E = cellfun(@(x) x.Control.E, states); 
I = cellfun(@(x) x.Control.I, states);
Tmax = cellfun(@(x) max(x.ThermalModel.T), states);
% [SOCN, SOCP] =  cellfun(@(x) model.calculateSOC(x), states);
time = cellfun(@(x) x.time, states); 

plot(time, E);

%%  energy density for discharge phase

E = cellfun(@(x) x.Control.E, dischargeStates); 
I = cellfun(@(x) x.Control.I, dischargeStates);
t = cellfun(@(x) x.time, dischargeStates);
mass = computeCellMass(model);
vol = sum(model.G.cells.volumes);
[Emid, Imid, energyDensity, specificEnergy, energy] = computeEnergyDensity(E, I, t, vol, mass);

figure
plot(energyDensity, Emid);
xlabel('Energy Density [Wh/L]');
ylabel('Voltage [V]');

figure
plot(specificEnergy, Emid);
xlabel('Specific Energy [Wh/kg]');
ylabel('Voltage [V]');
