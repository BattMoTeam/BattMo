
% clear the workspace and close open figures
clear all
close all

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa

%% shortcuts

ne = 'NegativeElectrode';
pe = 'PositiveElectrode';

am = 'ActiveMaterial';
gr = 'Graphite';
si = 'Silicon';

sd  = 'SolidDiffusion';
itf = 'Interface';

%% Setup the properties of Li-ion battery materials and cell design
jsonstruct = parseBattmoJson('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_silicon_graphite.json');

paramobj = SiliconGraphiteBatteryInputParams(jsonstruct);

rhoGr = paramobj.(ne).(am).(gr).(itf).density;
rhoSi = paramobj.(ne).(am).(si).(itf).density;

wfGr = 0.92; % weight fraction graphite
wfSi = 0.08; % weight fraction silicon

vfGr = wfGr/rhoGr;
vfSi = wfSi/rhoSi;
totV = (vfGr + vfSi);
vfGr = vfGr/totV;
vfSi = vfSi/totV;

paramobj.(ne).(am).(gr).activeMaterialFraction = vfGr;
paramobj.(ne).(am).(si).activeMaterialFraction = vfSi;

% paramobj.(ne).(am).(gr).(itf).theta0 = 0.01;
% paramobj.(ne).(am).(si).(itf).theta0 = 0.01;
% paramobj.(si).(sd).D0 = 1e-17;
% paramobj.(si).(itf).k0 = 1e-12;
% paramobj.(gr).(sd).D0 = 1e-16;

paramobj.scenario = 'first-charge';

paramobj = paramobj.validateInputParams();

gen = BatteryGenerator1D();

gen.xlength(4) = 1.8831*gen.xlength(4);

% gen.fac = 100;
% gen = gen.applyResolutionFactors();

% Now, we update the paramobj with the properties of the mesh. 
paramobj = gen.updateBatteryInputParams(paramobj);

model = SiliconGraphiteBattery(paramobj);
model.AutoDiffBackend= AutoDiffBackend();

inspectgraph = false;
if inspectgraph
    cgt = ComputationalGraphTool(model);
    [g, edgelabels] = cgt.getComputationalGraph('includeNodeNames', '[^d]T$', 'oneParentOnly', true);
    figure
    % h = plot(g, 'edgelabel', edgelabels, 'nodefontsize', 10);
    h = plot(g, 'nodefontsize', 10);
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
    switch model.scenario
      case 'discharge'
        inputI = model.Control.Imax;
        inputE = model.Control.lowerCutoffVoltage;
      case {'charge', 'first-charge'}
        inputI = -model.Control.Imax;
        inputE = model.Control.upperCutoffVoltage;
      otherwise
        error('initCase not recognized')
    end
    srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, inputI, inputE);
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

model.verbose = true;

nls = NonLinearSolver;
nls.errorOnFailure = false;

[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 


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
