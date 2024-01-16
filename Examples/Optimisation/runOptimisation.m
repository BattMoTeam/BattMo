clear all
close all

mrstModule add optimization

ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
am      = 'ActiveMaterial';
elyte   = 'Electrolyte';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
sep     = 'Separator';

% jsonstruct_geometry       = parseBattmoJson(fullfile(getParamsDir(), 'optim_1D_geometry.json'));
% jsonstruct_thicknesses    = parseBattmoJson(fullfile(getParamsDir(), 'optim_thicknesses.json'));
% jsonstruct_discretization = parseBattmoJson(fullfile(getParamsDir(), 'optim_1D_discretization.json'));
% jsonstruct_material       = parseBattmoJson(fullfile(getParamsDir(), 'optim_cell.json'));
% jsonstruct_timestepping   = parseBattmoJson(fullfile(getParamsDir(), 'optim_timestepping.json'));

% % we changed k0 to the value given by Sturm
% jsontruct_extra.(ne).(am).(itf).k0 = 2.22e-11;

% jsonstruct = mergeJsonStructs({jsontruct_extra          , ...
%                                jsonstruct_geometry      , ...
%                                jsonstruct_thicknesses   , ...
%                                jsonstruct_discretization, ...
%                                jsonstruct_material      , ...
%                                jsonstruct_timestepping});
% jsonstruct.Control.CRate = 1;

jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
jsonstruct.include_current_collectors = false;
jsonstruct.use_thermal = false;

jsonfilename = fullfile('Examples', 'JsonDataFiles', 'geometry1d.json');
jsonstruct_geometry = parseBattmoJson(jsonfilename);
jsonstruct = mergeJsonStructs({jsonstruct, jsonstruct_geometry});

output = runBatteryJson(jsonstruct, 'runSimulation', false, 'includeGridGenerator', true);

model   = output.model;
gridgen = output.gridGenerator;

css = CellSpecificationSummary(model);
css = css.updatePackingMass(60*gram);

%% setup default schedule
% We chose a time that should cover the battery life for all the parameters that will be tried.

CRate = model.Control.CRate;

totalTime = 1.4*hour/CRate;

% rampup stage
n  = 10;
dt = [];
dt = [dt; repmat(1e-4, n, 1).*1.5.^(1 : n)'];

% discharge stage
n     = 100;
dt    = [dt; repmat(totalTime/n, n, 1)];

times = [0; cumsum(dt)];
dt    = diff(times);
step  = struct('val', dt, 'control', ones(numel(dt), 1));

tup = 0.1; % rampup value for the current function, see rampupSwitchControl
srcfunc = @(time, I, E, Imax, lowerCutoffVoltage) rampupSwitchControl(time, ....
                                                                      tup , ...
                                                                      I   , ...
                                                                      E   , ...
                                                                      Imax, ...
                                                                      lowerCutoffVoltage);

control = struct('src', srcfunc, 'IEswitch', true);


schedule = struct('control', control, 'step', step);

%% Setup the initial state of the model

state0 = model.setupInitialState();

%% Setup the properties of the nonlinear solver

nls = NonLinearSolver();
nls.maxIterations  = 10;
nls.errorOnFailure = false;

model.nonlinearTolerance = 1e-5*model.Control.Imax;
model.verbose            = true;

dorunstartup = false;

if dorunstartup

    [wellSols, states, report] = simulateScheduleAD(state0, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);


    ind = cellfun(@(x) not(isempty(x)), states);
    states = states(ind);
    E = cellfun(@(x) x.(ctrl).E, states);
    I = cellfun(@(x) x.(ctrl).I, states);
    time = cellfun(@(x) x.time, states);

    figure
    plot(time, E)

    return

end

%%

simulatorSetup = struct('model', model, 'schedule', schedule, 'state0', state0);

parameters = {};

paramsetter = NPlengthSetter1D(model, gridgen, 1.1);

getValues = @(model, notused) paramsetter.getValues(model);
setValues = @(model, notused, v) paramsetter.setValues(model, v);

simulatorSetup.model = setValues(simulatorSetup.model, [], [80e-6; 0.3; 0.3]);
boxLims = [[50, 110]*1e-6; ...
           [0.15, 0.4]  ; ...
           [0.15, 0.4]];

css = css.updateModel(simulatorSetup.model);
css.printSpecifications();

% simulatorSetup.model = setValues(simulatorSetup.model, [], [75e-6; 0.15; 0.15]);

parameters{end+1} = ModelParameter(simulatorSetup            , ...
                                   'name'     , 'lengthPoros', ...
                                   'belongsTo', 'model'      , ...
                                   'location' , {''}         , ...
                                   'boxLims'  , boxLims      , ...
                                   'getfun'   , getValues    , ...
                                   'setfun'   , setValues);

params.E0    = 3.65;
params.alpha = 100;

mass = computeCellMass(model);
params.extraMass = css.packingMass;

objmatch = @(model, states, schedule, varargin) SpecificEnergyOutput(model, states, schedule, params, varargin{:});

options = {'NonLinearSolver', nls, 'OutputMinisteps', false};

% setup objective function scaling. We use initial value
p = getScaledParameterVector(simulatorSetup, parameters);
v0 = evalObjectiveBattmo(p, objmatch, simulatorSetup, parameters, 'GradientMethod', 'None', options{:});

options = horzcat(options, {'objScaling', v0});

doOptimization = true;

if doOptimization

    p_base = getScaledParameterVector(simulatorSetup, parameters);
    obj = @(p) evalObjectiveBattmo(p, objmatch, simulatorSetup, parameters, 'GradientMethod', 'AdjointAD', options{:});
    [v, p_opt, history] = unitBoxBFGS(p_base, obj, 'gradTol', 1e-4, 'objChangeTol', 1e-6);

    optimSetup = updateSetupFromScaledParameters(simulatorSetup, parameters, p_opt);

    vals = parameters{1}.getParameter(optimSetup);

end

doCompareGradient = true;

if doCompareGradient

    p = getScaledParameterVector(simulatorSetup, parameters);
    [vad, gad]   = evalObjectiveBattmo(p, objmatch, simulatorSetup, parameters, 'GradientMethod', 'AdjointAD', options{:});
    perturbationSize = {1e-10, 1e-5};
    [vnum, gnum] = evalObjectiveBattmo(p, objmatch, simulatorSetup, parameters, ...
                                       'GradientMethod', 'PerturbationADNUM'  , ...
                                       'PerturbationSize', perturbationSize   , ...
                                       options{:});

    fprintf('Gradient computed using adjoint:\n');
    display(gad);
    fprintf('Numerical gradient:\n');
    display(gnum);

end
