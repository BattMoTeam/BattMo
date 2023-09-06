%% Pseudo-Two-Dimensional (P2D) Lithium-Ion Battery Model
% This example demonstrates how to setup a P2D model of a Li-ion battery
% and run a simple simulation.

% clear the workspace and close open figures
clear all
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

gridgen.faceArea = 1;

% Now, we update the paramobj with the properties of the mesh. 
[paramobj, gridgen] = gridgen.updateBatteryInputParams(paramobj);

%%  Initialize the battery model. 

model = Battery(paramobj);

css = CellSpecificationSummary(model);
css.printSpecifications();

model.AutoDiffBackend= AutoDiffBackend();

%% Compute the nominal cell capacity and choose a C-Rate

CRate = model.Control.CRate;

total = 1.2*hour/CRate;

n    = 60;
dt   = total*1.4/n;
step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

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
% nls.timeStepSelector=StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);
model.nonlinearTolerance = 1e-5*model.Control.Imax;
% Set verbosity
model.verbose = true;


%% Set optimization parameters


dosometest = false;

if dosometest

    lsr = LengthSetter1D(gridgen, {ne, pe});
    
    reflengths = lsr.reflengths;

    modelAD = model;
    
    v = reflengths([1; 3]);
    v = initVariablesADI(v);

    modelAD = lsr.setLengths(modelAD, v);

    v = lsr.getLengths(model)
    vAD = lsr.getLengths(modelAD)
    
    mass = computeCellMass(model);
    cap  = computeCellCapacity(model);

    massAD = computeCellMass(modelAD);
    capAD  = computeCellCapacity(modelAD);

    return
    
end

dosometest = false;

if dosometest


    psr = PorositySetter(model, {ne, pe});
    
    refvalues = psr.refvalues;
    v = psr.getPorosities(model);

    v = initVariablesADI(v);

    modelAD = model;
    modelAD = psr.setPorosities(modelAD, v);

    v = psr.getPorosities(modelAD);
    
    return
end


% ind = cellfun(@(x) not(isempty(x)), states); 
% states = states(ind);

% E    = cellfun(@(x) x.Control.E, states); 
% time = cellfun(@(x) x.time, states); 

% plot(time, E, '*-');

%%

state0 = initstate;
simulatorSetup = struct('model', model, 'schedule', schedule, 'state0', state0);

%%

parameters= {};

includeLength = true;

if includeLength

    % lengthsetter = LengthSetter1D(gridgen, {ne, pe});
    lengthsetter = NPlengthSetter1D(model, gridgen, 1.1);
    
    getlength = @(model, dummy) lengthsetter.getLengths(model);
    setlength = @(model, dummy, v) lengthsetter.setLengths(model, v);

    parameters{end+1} = ModelParameter(simulatorSetup, ...
                                       'name'     , 'length'     , ...
                                       'belongsTo', 'model'      , ...
                                       'location' , {''}         , ...
                                       'boxLims'  , [40 100]*1e-6, ...
                                       'getfun'   , getlength    , ...
                                       'setfun'   , setlength);
end

includePorosity = true;

if includePorosity

    porosetter = PorositySetter(model, {ne, pe});
    
    getporo = @(model, dummy) porosetter.getPorosities(model);
    setporo = @(model, dummy, v) porosetter.setPorosities(model, v);
    
    parameters{end+1} = ModelParameter(simulatorSetup, ...
                                       'name'     , 'porosities', ...
                                       'belongsTo', 'model'     , ...
                                       'location' , {''}        , ...
                                       'boxLims'  , [0.1 0.9]   , ...
                                       'getfun'   , getporo     , ...
                                       'setfun'   , setporo);
    
end

% setup params

% params for cut-off function
params.E0    = 3;
params.alpha = 100;

%
mass = computeCellMass(model);
params.extraMass = 0.1*mass;

objmatch = @(model, states, schedule, varargin) SpecificEnergyOutput(model, states, schedule, params, varargin{:});
% objmatch = @(model, states, schedule, varargin) EnergyOutput(model, states, schedule, varargin{:});

options = {'NonLinearSolver', nls, 'OutputMinisteps', false};

% setup objective function scaling. We use initial value
p = getScaledParameterVector(simulatorSetup, parameters);
v0 = evalObjectiveBattmo(p, objmatch, simulatorSetup, parameters, 'GradientMethod', 'None', options{:});

options = horzcat(options, {'objScaling', v0});

dosmalltest = false;
if dosmalltest
    val = objmatch(model, states, schedule);
    sum([val{:}])
    v = parameters{1}.getParameterValue(simulatorSetup)
    obj = @(p) evalObjectiveBattmo(p, objmatch, simulatorSetup, parameters, options{:});
    vs = parameters{1}.scale(v)
    b = obj(vs)
end


doOptimization = true;

if doOptimization
    
    p_base = getScaledParameterVector(simulatorSetup, parameters);
    obj = @(p) evalObjectiveBattmo(p, objmatch, simulatorSetup, parameters, 'GradientMethod', 'AdjointAD', options{:});
    [v, p_opt, history] = unitBoxBFGS(p_base, obj, 'gradTol', 1e-4, 'objChangeTol', 1e-4);

end

doCompareGradient = false;

if doCompareGradient
    
    p = getScaledParameterVector(simulatorSetup, parameters);
    [vad, gad]   = evalObjectiveBattmo(p, objmatch, simulatorSetup, parameters, 'GradientMethod', 'AdjointAD', options{:});
    perturbationSize = {};
    if includeLength
        perturbationSize = horzcat(perturbationSize, {1e-10});
    end
    if includePorosity
        perturbationSize = horzcat(perturbationSize, {1e-5});
    end
    [vnum, gnum] = evalObjectiveBattmo(p, objmatch, simulatorSetup, parameters, ...
                                       'GradientMethod', 'PerturbationADNUM'  , ...
                                       'PerturbationSize', perturbationSize   , ...
                                       options{:});

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
