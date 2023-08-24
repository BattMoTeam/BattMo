%% Pseudo-Two-Dimensional (P2D) Lithium-Ion Battery Model
% This example demonstrates how to setup a P2D model of a Li-ion battery
% and run a simple simulation.

% clear the workspace and close open figures
% clear
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

paramobj = BatteryInputParams(jsonstruct);

%% Setup the geometry and computational mesh

gridgen = BatteryGenerator1D();

% Now, we update the paramobj with the properties of the mesh. 
[paramobj, gridgen] = gridgen.updateBatteryInputParams(paramobj);

%%  Initialize the battery model. 

model = Battery(paramobj);
model.AutoDiffBackend= AutoDiffBackend();

%% Compute the nominal cell capacity and choose a C-Rate

CRate = model.Control.CRate;

total = 1.2*hour/CRate;

n     = 40;
dt    = total*1.4/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

tup = 0.1; % rampup value for the current function, see rampupSwitchControl
srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                            model.Control.Imax, ...
                                            model.Control.lowerCutoffVoltage);

control = struct('src', srcfunc, 'IEswitch', true);


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
nls.timeStepSelector=StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);
model.nonlinearTolerance = 1e-3*model.Control.Imax;
% Set verbosity
model.verbose = true;

%% Run the simulation

% [~, states, ~] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

%% Process output and recover the output voltage and current from the output states.

lsr = LengthSetter1D(gridgen, {ne, pe});

dosometest = false;

if dosometest
    reflengths = lsr.reflengths;
    v = reflengths([1; 3]);
    v = initVariablesADI(v);

    model = lsr.setLength(model, v);
    v = lsr.getLength(model);

    mass = computeCellMass(model);
    cap = computeCellCapacity(model);
end

% ind = cellfun(@(x) not(isempty(x)), states); 
% states = states(ind);

% E    = cellfun(@(x) x.Control.E, states); 
% time = cellfun(@(x) x.time, states); 

% plot(time, E, '*-');

%%

state0 = initstate;
SimulatorSetup = struct('model', model, 'schedule', schedule, 'state0', state0);

%%

lengthsetter = LengthSetter1D(gridgen, {ne, pe});

getlength = @(model, dummy) lengthsetter.getLength(model);
setlength = @(model, dummy, v) lengthsetter.setLength(model, v);

parameters= {};
parameters{end+1} = ModelParameter(SimulatorSetup, ...
                                   'name'     , 'length'     , ...
                                   'belongsTo', 'model'      , ...
                                   'location' , {''}      , ...
                                   'boxLims'  , [40 100]*1e-6, ...
                                   'getfun'   , getlength   , ...
                                   'setfun'   , setlength);


cutoffVoltage = 3;
objmatch = @(model, states, schedule, varargin) SpecificEnergyOutput(model, states, schedule, cutoffVoltage, varargin{:});
a = objmatch(model, states, schedule);
fn       = @plotAfterStepIV;
obj      = @(p) evalObjectiveBattmo(p, objmatch, SimulatorSetup, parameters, 'objScaling', totval, 'afterStepFn', fn, 'NonLinearSolver', nls, 'OutputMinisteps', true);

v = parameters{1}.getParameterValue(SimulatorSetup)
SimulatorSetup = parameters{1}.setParameterValue(SimulatorSetup, v)
v = parameters{1}.getParameterValue(SimulatorSetup)
vs = parameters{1}.scale(v)
v = parameters{1}.unscale(vs)
b = evalObjectiveBattmo(vs, objmatch, SimulatorSetup, parameters, 'NonLinearSolver', nls, 'OutputMinisteps', true)
v = parameters{1}.getParameterValue(SimulatorSetup)

return

%%

doOptimization = true;

if doOptimization
    
    p_base = getScaledParameterVector(SimulatorSetup, parameters);
    p_base = p_base - 0.1;
    [v, p_opt, history] = unitBoxBFGS(p_base, obj, 'gradTol', 1e-7, 'objChangeTol', 1e-4);
    
end


doCompareGradient = false;
if doCompareGradient
    
    p = getScaledParameterVector(SimulatorSetup, parameters);
    [vad, gad]   = evalObjectiveBattmo(p, objmatch, SimulatorSetup, parameters, 'Gradient', 'AdjointAD');
    [vnum, gnum] = evalObjectiveBattmo(p, objmatch, SimulatorSetup, parameters, 'Gradient', 'PerturbationADNUM', 'PerturbationSize', 1e-5);

    fprintf('Gradient computed using adjoint:\n');
    display(gad);
    fprintf('Numerical gradient:\n');
    display(gnum);
    
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
