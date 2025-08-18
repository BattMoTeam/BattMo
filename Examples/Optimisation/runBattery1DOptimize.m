%% Setup

% Clear the workspace and close open figures
clear all
close all

%% Setup the properties of Li-ion battery materials and cell design

jsonfilename = fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_nmc_graphite.json');
jsonstruct_material = parseBattmoJson(jsonfilename);

jsonfilename = fullfile('Examples', 'JsonDataFiles', 'geometry1d.json');
jsonstruct_geometry = parseBattmoJson(jsonfilename);

jsonfilename = fullfile('Examples', 'JsonDataFiles', 'cc_discharge_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);

jsonstruct = mergeJsonStructs({jsonstruct_geometry , ...
                               jsonstruct_material , ...
                               jsonstruct_control});

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

jsonstruct.(ne).(am).diffusionModelType = 'simple';
jsonstruct.(pe).(am).diffusionModelType = 'simple';

jsonstruct.(ctrl).useCVswitch = true;

output = runBatteryJson(jsonstruct, 'runSimulation', false);

simsetup = output.simsetup;

%% Setup the nonlinear solver

nls = NonLinearSolver();

% Change the number of maximum nonlinear iterations
nls.maxIterations = 10;

% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false;

simsetup.NonLinearSolver = nls;
% Change tolerance for the nonlinear iterations
simsetup.model.nonlinearTolerance = 1e-3*simsetup.model.Control.Imax;

% Set verbosity
simstup.model.verbose = false;

states = simsetup.run;

%% Process output and recover the output voltage and current from the output states.

ind = cellfun(@(x) not(isempty(x)), states);
states = states(ind);

E    = cellfun(@(x) x.Control.E, states);
I    = cellfun(@(x) x.Control.I, states);
time = cellfun(@(x) x.time, states);

doPlot = false;

if doPlot
    figure;
    plot(time/hour, E, '*-', 'displayname', 'initial');
    xlabel('time  / h');
    ylabel('voltage  / V');
    grid on
end

%% Calculate the energy

model    = simsetup.model;
schedule = simsetup.schedule;

obj = @(simsetup, states, varargin) EnergyOutput(simsetup, states, varargin{:});
vals = obj(simsetup, states);
totval = sum([vals{:}]);

% Compare with trapezoidal integral: they should be about the same
totval_trapz = trapz(time, E.*I);
fprintf('Rectangle rule: %g Wh, trapezoidal rule: %g Wh\n', totval/hour, totval_trapz/hour);

%% Setup the optimization problem

parameters = {};

paramsetter = PorositySetter(model, {ne, sep, pe});

getporo = @(model, notused) paramsetter.getValues(model);
setporo = @(model, notused, v) paramsetter.setValues(model, v);

parameters = addParameter(parameters , simsetup  , ...
                          'name'     , 'porosity', ...
                          'belongsTo', 'model'   , ...
                          'boxLims'  , [0.1      , 0.9], ...
                          'location' , {''}      , ...
                          'getfun'   , getporo   , ...
                          'setfun'   , setporo);

getfun = @(model, notused) model.Control.Imax;
setfun = @(model, notused, v) setImax(model, v);

parameters = addParameter(parameters , simsetup               , ...
                          'name'     , 'Imax'                 , ...
                          'belongsTo', 'model'                , ...
                          'boxLims'  , model.Control.Imax*[0.5, 2], ...
                          'location' , {''}                   , ...
                          'getfun'   , getfun                 , ...
                          'setfun'   , setfun);

%% Setup the objective function and auxiliary plotting

objmatch = @(simsetup, states, varargin) EnergyOutput(simsetup, states, varargin{:});
if doPlot
    fn = @plotAfterStepIV;
else
    fn = [];
end
obj = @(p) evalObjectiveBattmo(p, objmatch, simsetup, parameters, 'objScaling', totval, 'afterStepFn', fn);

%% Setup initial parameters

% The parameters must be scaled to [0,1]
p_base = getScaledParameterVector(simsetup, parameters);
p_base = max(0, p_base - 0.1);

%% Optimize

% Solve the optimization problem using BFGS. One can adjust the
% tolerances and the maxIt option to see how it effects the
% optimum.
[v, p_opt, history] = unitBoxBFGS(p_base, obj, 'gradTol', 1e-7, 'objChangeTol', 1e-4, 'maxIt', 20);

return
% Compute objective at optimum
opt_simsetup = updateSetupFromScaledParameters(simsetup, parameters, p_opt);
[~, states_opt, ~] = simulateScheduleAD(opt_simsetup.state0, opt_simsetup.model, opt_simsetup.schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);
time_opt = cellfun(@(x) x.time, states_opt);
E_opt = cellfun(@(x) x.Control.E, states_opt);
I_opt = cellfun(@(x) x.Control.I, states_opt);
totval_trapz_opt = trapz(time_opt, E_opt.*I_opt);

% Print optimal parameters
fprintf('Base and optimized parameters:\n');
for k = 1:numel(parameters)
    % Get the original and optimized values
    p0 = parameters{k}.getParameter(simsetup);
    pu = parameters{k}.getParameter(opt_simsetup);

    % Print
    fprintf('%s\n', parameters{k}.name);
    fprintf('%g %g\n', p0, pu);
end

fprintf('Energy changed from %g to %g mWh\n', totval_trapz/hour/milli, totval_trapz_opt/hour/milli);

%% Plot

if doPlot
    % Plot
    figure; hold on; grid on
    E    = cellfun(@(x) x.Control.E, states);
    time = cellfun(@(x) x.time, states);
    plot(time/hour, E, '*-', 'displayname', 'initial');
    plot(time_opt/hour, E_opt, 'r*-', 'displayname', 'optimized');
    xlabel('time  / h');
    ylabel('voltage  / V');
    legend;
end

%%
doCompareGradient = false;
if doCompareGradient

    p = getScaledParameterVector(simsetup, parameters);
    [vad, gad]   = evalObjectiveBattmo(p, objmatch, simsetup, parameters, 'gradientMethod', 'AdjointAD');
    [vnum, gnum] = evalObjectiveBattmo(p, objmatch, simsetup, parameters, 'gradientMethod', 'PerturbationADNUM', 'PerturbationSize', 1e-5);

    fprintf('Gradient computed using adjoint:\n');
    display(gad);
    fprintf('Numerical gradient:\n');
    display(gnum);
    fprintf('Relative error:\n')
    display(abs(gnum-gad)./abs(gad));
end

function model = setImax(model, Imax)

    model.Control.Imax = Imax;
    
end

%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
