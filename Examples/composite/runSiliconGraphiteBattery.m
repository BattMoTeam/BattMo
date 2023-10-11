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

% Now, we update the paramobj with the properties of the mesh. 
paramobj = gen.updateBatteryInputParams(paramobj);

% We instantiate the model
model = Battery(paramobj);

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

tup = 0.1; % rampup value for the current function, see rampupSwitchControl
srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                            model.Control.Imax, ...
                                            model.Control.lowerCutoffVoltage);
% we setup the control by assigning a source and stop function.
control = struct('src', srcfunc, 'IEswitch', true);

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

[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

dischargeStates = states;

%%

initstate = states{end};

srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                            -model.Control.Imax, ...
                                            model.Control.upperCutoffVoltage);
control = struct('src', srcfunc, 'IEswitch', true);
schedule = struct('control', control, 'step', step); 

[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

chargeStates = states;

%%

allStates = vertcat(dischargeStates, chargeStates); 

set(0, 'defaultlinelinewidth', 3);
set(0, 'DefaultAxesFontSize', 16);
set(0, 'defaulttextfontsize', 18);

E = cellfun(@(x) x.Control.E, allStates); 
I = cellfun(@(x) x.Control.I, allStates);
Tmax = cellfun(@(x) max(x.ThermalModel.T), allStates);
time = cellfun(@(x) x.time, allStates); 

figure
plot(time/hour, E);
xlabel('Time / h');
ylabel('Voltage / V');
title('Voltage')

figure
plot(time/hour, I);
xlabel('Time / h');
ylabel('Current / I');
title('Current')

figure
hold on

for istate = 1 : numel(allStates)
    allStates{istate} = model.evalVarName(allStates{istate}, {ne, co, 'SOC'});
end

SOC  = cellfun(@(x) x.(ne).(co).SOC, allStates); 
SOC1 = cellfun(@(x) x.(ne).(co).(am1).SOC, allStates);
SOC2 = cellfun(@(x) x.(ne).(co).(am2).SOC, allStates);

plot(time/hour, SOC, 'displayname', 'SOC - cumulated');
plot(time/hour, SOC1, 'displayname', 'SOC - Graphite');
plot(time/hour, SOC2, 'displayname', 'SOC - Silicon');

xlabel('Time / h');
ylabel('SOC / -');
title('SOCs')

legend show
