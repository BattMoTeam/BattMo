%% Setup
mrstModule add ad-core optimization mpfa mrst-gui

clear all
close all

ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
am      = 'ActiveMaterial';
co      = 'Coating';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
sep     = 'Separator';

%% Choose battery type

jsonParams  = parseBattmoJson(fullfile('ParameterData'        , ...
                                       'BatteryCellParameters', ...
                                       'LithiumIonBatteryCell', ...
                                       'lithium_ion_battery_nmc_graphite.json'));
jsonGeom    = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometry1d.json'));
jsonControl = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'cc_discharge_control.json'));
jsonSim     = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'simulation_parameters.json'));

json = mergeJsonStructs({jsonParams, jsonGeom, jsonControl, jsonSim});

json.use_thermal = false;

json.Control.useCVswitch = true;

% % Test finer time discretization
% json.TimeStepping.numberOfTimeSteps = 80;
% json.TimeStepping.numberOfRampupSteps = 10;

% Optionally validate the json struct
validateJson = false;

%% Run with initial guess
json0 = json;
output0 = runBatteryJson(json0, 'validateJson', validateJson);

simSetup = struct('model'   , output0.model   , ...
                  'schedule', output0.schedule, ...
                  'state0'  , output0.initstate);

%% Setup parameters to be optimized
params = [];

% Electrolyte bryggeman coefficient

params = addParameter(params                   , ...
                      simSetup                 , ...
                      'name', 'elyte_bruggeman', ...
                      'belongsTo', 'model'     , ...
                      'boxLims', [1, 3]        , ...
                      'location', {elyte, 'bruggemanCoefficient'});

% Exchange current densities in the Butler-Volmer eqn

params = addParameter(params                     , ...
                      simSetup                   , ...
                      'name'     , 'ne_k0'       , ...
                      'belongsTo', 'model'       , ...
                      'boxLims'  , [1e-12 , 1e-10], ...
                      'location' , {ne, co, am, itf, 'reactionRateConstant'});

params = addParameter(params                  , ...
                      simSetup                , ...
                      'name', 'pe_k0'         , ...
                      'belongsTo', 'model'    , ...
                      'boxLims', [1e-12, 1e-10], ...
                      'location', {pe, co, am, itf, 'reactionRateConstant'});

% Volumetric surface areas

params = addParameter(params               , ...
                      simSetup             , ...
                      'name', 'ne_vsa'     , ...
                      'belongsTo', 'model' , ...
                      'boxLims', [1e5, 1e6], ...
                      'location', {ne, co, am, itf, 'volumetricSurfaceArea'});

params = addParameter(params               , ...
                      simSetup             , ...
                      'name', 'pe_vsa'     , ...
                      'belongsTo', 'model' , ...
                      'boxLims', [1e5, 1e6], ...
                      'location', {pe, co, am, itf, 'volumetricSurfaceArea'});

%% Generate "experimental" data that we want to match

expSimSetup = simSetup;

for ip = 1:numel(params)
    loc         = params{ip}.location;
    new         = mean(params{ip}.boxLims);
    expSimSetup = params{ip}.setParameter(expSimSetup, new);
end

%%

simulator = Simulator(expSimSetup.model                 , ...
                      'initialState', expSimSetup.state0, ...
                      'schedule', expSimSetup.schedule);

expStates = simulator.run();

%% Setup the optimization problem

% Objective function
objective = @(model, states, schedule, varargin) leastSquaresEI(expSimSetup.model, states, expStates, expSimSetup.schedule, varargin{:});

% Debug: the objective function evaluated at the experimental values
% should be zero
objval = objective(expSimSetup.model,expStates, expSimSetup.schedule);
assert(max(abs([objval{:}])) == 0.0);

% Function for gradient computation
objVals           = objective(output0.model, output0.states, output0.schedule);
objScaling        = sum([objVals{:}]);
objectiveGradient = @(p) evalObjectiveBattmo(p, objective, simSetup, params, 'objScaling', objScaling);

%% Optional debug: Compare gradients using adjoints and finite difference approximation
debug = true
if debug
    pTmp = getScaledParameterVector(simSetup, params);

    [vad, gad, tmp, statesx, setupNewx, lambdas] = evalObjectiveBattmo(pTmp, objective, simSetup, params, ...
                                     'gradientMethod', 'AdjointAD');

    [vnum, gnum] = evalObjectiveBattmo(pTmp, objective, simSetup, params, ...
                                       'gradientMethod', 'PerturbationADNUM');

    fprintf('Adjoint and finite difference derivatives and the relative error:\n');
    disp([gad, gnum, abs(gad-gnum)./abs(gad)])

    plotDashboardAdjoint(setupNewx.model, lambdas, 'step', 1);
    plotDashboardAdjoint(setupNewx.model, lambdas, 'step', 0);

    return
end

%% Optimize

p0scaled     = getScaledParameterVector(simSetup, params);
gradTol      = 1e-7;
objChangeTol = 1e-6;
maxIt        = 25;

[v, pOptTmp, history] = unitBoxBFGS(p0scaled                    , ....
                                    objectiveGradient           , ...
                                    'maximize'    , false       , ...
                                    'gradTol'     , gradTol     , ...
                                    'objChangeTol', objChangeTol, ...
                                    'maxIt'       , maxIt       , ...
                                    'logplot'     , true);
numIt = numel(history.val);

%% Unscale

optimSetup = simSetup;

pOpt = zeros(size(pOptTmp));
for k = 1 : numel(params)
    pOpt(k)    = params{k}.unscale(pOptTmp(k));
    optimSetup = params{k}.setParameter(optimSetup, pOpt(k));
end


%% Run model with pOpt params

simulator = Simulator(optimSetup.model, ...
                      'initialState', optimSetup.state0, ...
                      'schedule', optimSetup.schedule);

optimStates = simulator.run();

%% plotting

do_plot = true;
if do_plot
    set(0, 'defaultlinelinewidth', 2)

    getTime = @(states) cellfun(@(state) state.time, states);
    getE = @(states) cellfun(@(state) state.Control.E, states);

    t0   = getTime(output0.states);
    E0   = getE(output0.states);
    tOpt = getTime(optimStates);
    EOpt = getE(optimStates);
    tExp = getTime(expStates);
    EExp = getE(expStates);

    h = figure; hold on; grid on; axis tight
    plot(t0/hour, E0, 'displayname', 'E_{0}')
    plot(tExp/hour, EExp, '--', 'displayname', 'E_{exp}');
    plot(tOpt/hour, EOpt, ':', 'displayname', 'E_{opt}')
    legend;

end

%% Summarize
pOrig = cellfun(@(p) p.getParameter(simSetup), params)';

fprintf('Initial guess:\n');
fprintf('%g\n', pOrig);

fprintf('Fitted values (* means we hit the box limit):\n');
tol = 1e-3;
for k = 1:numel(params)
    hit = '';
    if abs(pOptTmp(k)) < tol || abs(pOptTmp(k)-1) < tol
        hit = '*';
    end
    fprintf('%g %s\n', pOpt(k), hit);
end

pExp = cellfun(@(p) p.getParameter(expSimSetup), params)';

fprintf('\nExperimental values:\n');
fprintf('%g\n', pExp);

fprintf('\nIterations:\n')
fprintf('%g\n', numIt);


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
