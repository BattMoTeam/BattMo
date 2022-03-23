%% Pseudo-Two-Dimensional (P2D) Lithium-Ion Battery Model
% This example demonstrates how to setup a P2D model of a Li-ion battery
% and run a simple simulation.

% clear the workspace and close open figures
clear
close all
clc


%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa

%% Setup the properties of Li-ion battery materials and cell design
% The properties and parameters of the battery cell, including the
% architecture and materials, are set using an instance of
% :class:`BatteryInputParams <Battery.BatteryInputParams>`. This class is
% used to initialize the simulation and it propagates all the parameters
% throughout the submodels. The input parameters can be set manually or
% provided in json format. All the parameters for the model are stored in
% the paramobj object.
jsonstruct = parseBattmoJson('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json');
% jsonstruct.Control.controlPolicy = 'CCCV';
paramobj = BatteryInputParams(jsonstruct);
% paramobj.SOC = 0.02;

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

%% Setup the geometry and computational mesh
% Here, we setup the 1D computational mesh that will be used for the
% simulation. The required discretization parameters are already included
% in the class BatteryGenerator1D. 
gen = BatteryGenerator1D();
gen.fac = 1;
gen = gen.applyResolutionFactors();

% Now, we update the paramobj with the properties of the mesh. 
paramobj = gen.updateBatteryInputParams(paramobj);

% !!! REMOVE THIS. SET THE RIGHT VALUES IN THE JSON !!! In this case, we
% change some of the values of the paramaters that were given in the json
% file to other values. This is done directly on the object paramobj. 
paramobj.(ne).(cc).EffectiveElectricalConductivity = 1e5;
paramobj.(pe).(cc).EffectiveElectricalConductivity = 1e5;
paramobj.(thermal).externalTemperature = paramobj.initT;

%%  Initialize the battery model. 
% The battery model is initialized by sending paramobj to the Battery class
% constructor. see :class:`Battery <Battery.Battery>`.
model = Battery(paramobj);
model.AutoDiffBackend= AutoDiffBackend();

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
    total = 1*hour/CRate;
  otherwise
    error('control policy not recognized');
end

n     = 20;
dt    = total/n;
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
%%
nc = 3;
nst = numel(step.control)
ind = floor(([0:nst-1]/nst)*nc)+1
%%
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
NLS.errorOnFailure = false; 
%nls.timeStepSelector=StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);
% Change default tolerance for nonlinear solver
nls.timeStepSelector = SimpleTimeStepSelector();
model.nonlinearTolerance = 1e-3*model.Control.Imax;
% Set verbosity
model.verbose = true;

%% Run the simulation
[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

%% Process output and recover the output voltage and current from the output states.
ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
Enew = cellfun(@(x) x.Control.E, states); 
Inew = cellfun(@(x) x.Control.I, states);
Tmax = cellfun(@(x) max(x.ThermalModel.T), states);
[SOCN, SOCP] =  cellfun(@(x) model.calculateSOC(x), states);
time = cellfun(@(x) x.time, states); 
plot(time,Enew,'*-')
%% Plot the the output voltage and current
% gradients to control
obj = @(model, states, schedule, varargin) EnergyOutput(model, states, schedule,varargin{:});%,'step',step);
vals = obj(model, states, schedule)
totval = sum([vals{:}]);
%%
scaling = struct('boxLims',[],'obj',[]);
scaling.boxLims = model.Control.Imax*[0.9,1.1];
scaling.obj = 1;
%%
state0 = initstate;
%obj1 = @(step,state,model, states, schedule, varargin) EnergyOutput(model, states, schedule,varargin{:},'step',step);
f = @(u)evalObjectiveBattmo(u, obj, state0, model, schedule, scaling);
%schedule.control(:).Imax = model.Control.Imax;
u_base = battmoSchedule2control(schedule, scaling);
[val, gad] =evalObjectiveBattmo(u_base, obj, state0, model, schedule, scaling);
%%
[val,gnum] =evalObjectiveBattmo(u_base, obj, state0, model, schedule, scaling,'Gradient','numerical');
%%
return
% Get function handle for objective evaluation

%%
mrstModule add optimization
%% Run optimization with default options
% gradient to 
%schedule_opt = control2schedule(u_opt, schedule, scaling);
%model.toleranceCNV = 1e-6;
objsens =@(tstep,model, state) EnergyOutput(model, states, schedule,'tstep',tstep,'computePartials',true,'state',state);
%matchObservedOW(model, states, schedule, states_ref,...
%    'computePartials', true, 'tstep', tstep, weighting{:},'state',state);
SimulatorSetup = struct('model', model, 'schedule', schedule, 'state0', state0);
parameters = [];
property = 'Imax'
%getfun = @(x,location) x.Imax; 
setfun = @(x,location, v) struct('Imax',v,'src', @(time,I,E) rampupSwitchControl(time, 0.1, I, E, ...
                                                v, ...
                                                2.0),'IEswitch',true);
boxlims = model.Control.Imax*[0.9,1.1];
parameters={};
for i = 1:3
parameters{end+1} = ModelParameter(SimulatorSetup,'name','Imax',...
    'belongsTo','schedule',...
    'boxLims',boxlims,...
    'location',{'control','Imax'},...
    'getfun',[],'setfun',setfun,'controlSteps',i)
end

%%
%parameters = addParameter(parameters, SimulatorSetup,'belongeTo', 'name', property ,'type','multiplier');
raw_sens = computeSensitivitiesAdjointADBattmo(SimulatorSetup, states, parameters, objsens)
%%
ff =fieldnames(raw_sens)
gadnew = nan(numel(ff),1)
for i =1:numel(ff)
    gadnew(i) = raw_sens.(ff{i})
end
%%
gadnew*diff(boxlims)
%%
%objmatch = @(model, states, schedule, states_ref, compDer, tstep, state) ...
%    EnergyOutput(model, states, schedule,'tstep',tstep,'computePartials',true,'state',state);
objmatch = @(model, states, schedule, varargin) EnergyOutput(model, states, schedule,varargin{:});%
pvec = getScaledParameterVector(SimulatorSetup, parameters);
objh = @(p) evalMatchBattmo(p, objmatch, SimulatorSetup, parameters, []);
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 30 iterations
%%
pert = 0.1;
[vad,gad] = evalMatchBattmo(pvec+pert, objmatch, SimulatorSetup, parameters, 'Gradient','AdjointAD')
[vnum,gnum] = evalMatchBattmo(pvec+pert, objmatch, SimulatorSetup, parameters, 'Gradient','PerturbationADNUM','PerturbationSize',1e-4)
%%
%[v, p_opt, history] = unitBoxBFGS(pvec, objh, 'objChangeTol', 1e-5, ...
%    'maxIt', 30, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

%setup_opt = updateSetupFromScaledParameters(setup_init, parameters, p_opt); 
%[wellSols_opt, states_opt] = simulateScheduleAD(setup_opt.state0, setup_opt.model, setup_opt.schedule);

%% Run optimization with default options
% gradient to 
%schedule_opt = control2schedule(u_opt, schedule, scaling);
%model.toleranceCNV = 1e-6;
%objsens =@(tstep,model, state) EnergyOutput(model, states, schedule,'tstep',tstep,'computePartials',true,'state',state);
%matchObservedOW(model, states, schedule, states_ref,...
%    'computePartials', true, 'tstep', tstep, weighting{:},'state',state);
SimulatorSetup = struct('model', model, 'schedule', schedule, 'state0', state0);
parameters = [];
property = 'porevolume'
getfun = @(x,location) x.(location).volumeFraction; 
setfun = @(x,location, v) setfunctionWithName(x,location,v,'volumeFraction');
parameters={};
parameters{end+1} = ModelParameter(SimulatorSetup,'name','volumeFraction',...
    'belongsTo','model',...
    'location',{'Electrolyte','volumeFraction'},...
    'getfun',[],'setfun',[])

%parameters = addParameter(parameters, SimulatorSetup,'belongeTo', 'name', property ,'type','multiplier');
raw_sens = computeSensitivitiesAdjointADBattmo(SimulatorSetup, states, parameters, objsens);
%%
objmatch = @(model, states, schedule, varargin) EnergyOutput(model, states, schedule,varargin{:});%
pvec = getScaledParameterVector(SimulatorSetup, parameters);
objh = @(p) evalMatchBattmo(p, objmatch, SimulatorSetup, parameters, []);
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 30 iterations
%%
pert = 0.0;
[vad,gad] = evalMatchBattmo(pvec+pert, objmatch, SimulatorSetup, parameters, 'Gradient','AdjointAD')
[vnum,gnum] = evalMatchBattmo(pvec+pert, objmatch, SimulatorSetup, parameters, 'Gradient','PerturbationADNUM','PerturbationSize',1e-4)



%% 

%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
