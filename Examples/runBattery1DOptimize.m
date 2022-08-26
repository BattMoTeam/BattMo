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
paramobj.use_thermal = false;
paramobj.include_current_collectors = false;
paramobj.NegativeElectrode.ActiveMaterial.externalCouplingTerm = paramobj.NegativeElectrode.CurrentCollector.externalCouplingTerm;
paramobj.PositiveElectrode.ActiveMaterial.externalCouplingTerm = paramobj.PositiveElectrode.CurrentCollector.externalCouplingTerm;
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
%%
nc = 1;
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
% gradients to controlmodel.use_thermal =false;

obj = @(model, states, schedule, varargin) EnergyOutput(model, states, schedule,varargin{:});%,'step',step);
vals = obj(model, states, schedule)
totval = sum([vals{:}]);
if(false)
%%
scaling = struct('boxLims',[],'obj',[]);
scaling.boxLims = model.Control.Imax*[0.9,1.1];
scaling.obj = totval;
%%
%state0 = initstate;
%obj1 = @(step,state,model, states, schedule, varargin) EnergyOutput(model, states, schedule,varargin{:},'step',step);
f = @(u)evalObjectiveBattmo(u, obj, state0, model, schedule, scaling);
%schedule.control(:).Imax = model.Control.Imax;
u_base = battmoSchedule2control(schedule, scaling);
end
%[v, u_opt, history] = unitBoxBFGS(u_base, f);
%schedule_opt = control2schedule(u_opt, schedule, scaling);
%return
%%

% Get function handle for objective evaluation

%%
if(false)
    %%
mrstModule add optimization
SimulatorSetup = struct('model', model, 'schedule', schedule, 'state0', state0);
parameters = [];
property = 'Imax'
%getfun = @(x,location) x.Imax; 
setfun = @(x,location, v) struct('Imax',v,'src', @(time,I,E) rampupSwitchControl(time, 0.1, I, E, ...
                                                v, ...
                                                3.0),'IEswitch',true);
boxlims = model.Control.Imax*[0.5,1.5];
parameters={};
for i = 1:max(schedule.step.control)
parameters{end+1} = ModelParameter(SimulatorSetup,'name','Imax',...
    'belongsTo','schedule',...
    'boxLims',boxlims,...
    'location',{'control','Imax'},...
    'getfun',[],'setfun',setfun,'controlSteps',i)
end
fn =@ plotAfterStepIV
%
%objmatch = @(model, states, schedule, states_ref, compDer, tstep, state) ...
%    EnergyOutput(model, states, schedule,'tstep',tstep,'computePartials',true,'state',state);
objmatch = @(model, states, schedule, varargin) EnergyOutput(model, states, schedule,varargin{:});%
p_base = getScaledParameterVector(SimulatorSetup, parameters);
obj = @(p) evalMatchBattmo(p, objmatch, SimulatorSetup, parameters,'objScaling',-totval,'afterStepFn',fn);
%%
figure(33),clf
[v, p_opt, history] = unitBoxBFGS(p_base, obj);
%%
end
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 30 iterations
%%
state0 = initstate;
SimulatorSetup = struct('model', model, 'schedule', schedule, 'state0', state0);
parameters={};

part = ones(model.Electrolyte.Separator.G.cells.num,1)
parameters{end+1} = ModelParameter(SimulatorSetup,'name','porosity_sep',...
    'belongsTo','model',...
    'boxLims',[0.1 0.9],...
    'lumping',part,...
    'location',{'Electrolyte','Separator','porosity'},...
    'getfun',[],'setfun',[])
%parameters={};
part = ones(model.NegativeElectrode.ActiveMaterial.G.cells.num,1);
parameters{end+1} = ModelParameter(SimulatorSetup,'name','porosity_nam',...
    'belongsTo','model',...
    'boxLims',[0.07 0.5],...
    'lumping',part,...
    'location',{'NegativeElectrode','ActiveMaterial','porosity'},...
    'getfun',[],'setfun',[])
part = ones(model.NegativeElectrode.ActiveMaterial.G.cells.num,1);
parameters{end+1} = ModelParameter(SimulatorSetup,'name','porosity_pam',...
    'belongsTo','model',...
    'boxLims',[0.07 0.5],...
    'lumping',part,...
    'location',{'PositiveElectrode','ActiveMaterial','porosity'},...
    'getfun',[],'setfun',[])

setfun = @(x,location, v) struct('Imax',v,'src', @(time,I,E) rampupSwitchControl(time, 0.1, I, E, ...
                                                v, ...
                                                3.0),'IEswitch',true);
boxlims = model.Control.Imax*[0.5,1.5];
%parameters={};
for i = 1:max(schedule.step.control)
parameters{end+1} = ModelParameter(SimulatorSetup,'name','Imax',...
    'belongsTo','schedule',...
    'boxLims',boxlims,...
    'location',{'control','Imax'},...
    'getfun',[],'setfun',setfun,'controlSteps',i)
end
objmatch = @(model, states, schedule, varargin) EnergyOutput(model, states, schedule,varargin{:});%
p_base = getScaledParameterVector(SimulatorSetup, parameters);
fn =@ plotAfterStepIV

obj = @(p) evalMatchBattmo(p, objmatch, SimulatorSetup, parameters,'objScaling',-totval,'afterStepFn',fn);
figure(33),clf
%%
%[v, p_opt, history] = unitBoxBFGS(p_base-0.1, obj,'gradTol',1e-7,'objChangeTol',1e-4);

%% 
% check derivatives
return 
%%
% test adjoint
parameters_tmp = {parameters{3}}
parameters_tmp = parameters;%{parameters{3},parameters{2}};
p_tmp = getScaledParameterVector(SimulatorSetup, parameters_tmp);
pert = 0.0;
[vad,gad]   = evalMatchBattmo(p_tmp+pert, objmatch, SimulatorSetup, parameters_tmp, 'Gradient','AdjointAD')
[vnum,gnum] = evalMatchBattmo(p_tmp+pert, objmatch, SimulatorSetup, parameters_tmp, 'Gradient','PerturbationADNUM','PerturbationSize',1e-5)
gad
gnum
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