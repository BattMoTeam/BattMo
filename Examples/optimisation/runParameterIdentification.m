if mrstPlatform('octave')
    error('This demo cannot be run from Octave');
end

%% Setup
mrstModule add ad-core optimization mpfa mrst-gui

clear all
close all

do_plot = true;

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

jsonParams = parseBattmoJson(fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', 'lithium_ion_battery_nmc_graphite.json'));

% jsonParams = parseBattmoJson(fullfile('ParameterData', 'ParameterSets', 'Xu2015', 'lfp.json'));
% jsonGeom = parseBattmoJson(fullfile('Examples', 'jsondatafiles', 'geometry1d.json'));
% json = mergeJsonStructs({jsonParams, jsonGeom});

json = jsonParams;

%% Run with initial guess
json0 = json;
[states0, model, schedule, initstate] = runPImodel(json0);

simSetup = struct('model', model, ...
                  'schedule', schedule, ...
                  'state0', initstate);

%% Setup parameters to be optimized
params = [];

% Bruggeman coefficients
params = addParameter(params, simSetup, ...
                      'name', 'ne_co_bruggeman', ...
                      'belongsTo', 'model', ...
                      'boxLims', [1, 3], ...
                      'location', {ne, co, 'bruggemanCoefficient'});
params = addParameter(params, simSetup, ...
                      'name', 'elyte_bruggeman', ...
                      'belongsTo', 'model', ...
                      'boxLims', [1, 3], ...
                      'location', {elyte, 'bruggemanCoefficient'});
params = addParameter(params, simSetup, ...
                      'name', 'sep_bruggeman', ...
                      'belongsTo', 'model', ...
                      'boxLims', [1, 3], ...
                      'location', {sep, 'bruggemanCoefficient'});
params = addParameter(params, simSetup, ...
                      'name', 'pe_co_bruggeman', ...
                      'belongsTo', 'model', ...
                      'boxLims', [1, 3], ...
                      'location', {pe, co, 'bruggemanCoefficient'});

% Exchange current densities in the Butler-Volmer eqn
params = addParameter(params, simSetup, ...
                      'name', 'ne_k0', ...
                      'belongsTo', 'model', ...
                      'scaling', 'log', ...
                      'boxLims', [1e-12, 1e-9], ...
                      'location', {ne, co, am, itf, 'reactionRateConstant'});
params = addParameter(params, simSetup, ...
                      'name', 'pe_k0', ...
                      'belongsTo', 'model', ...
                      'scaling', 'log', ...
                      'boxLims', [1e-12, 1e-9], ...
                      'location', {pe, co, am, itf, 'reactionRateConstant'});

% Volumetric surface areas
params = addParameter(params, simSetup, ...
                      'name', 'ne_vsa', ...
                      'belongsTo', 'model', ...
                      'boxLims', [1e5, 1e7], ...
                      'location', {ne, co, am, itf, 'volumetricSurfaceArea'});
params = addParameter(params, simSetup, ...
                      'name', 'pe_vsa', ...
                      'belongsTo', 'model', ...
                      'boxLims', [1e5, 1e7], ...
                      'location', {pe, co, am, itf, 'volumetricSurfaceArea'});

%% Generate "experimental" data
jsonExp = json;
for ip = 1:numel(params)
    loc = params{ip}.location;
    orig = getfield(jsonExp, loc{:});
    new = mean(params{ip}.boxLims);
    jsonExp = setfield(jsonExp, loc{:}, new);
end
pExp = cellfun(@(p) getfield(jsonExp, p.location{:}), params)';
statesExp = runPImodel(jsonExp);

%% Setup the optimization problem

% Objective function
objective = @(model, states, schedule, varargin) leastSquaresEI(model, states, statesExp, schedule, varargin{:});

% Debug: the objective function evaluated at the experimental values should be zero
objval = objective(model, statesExp, schedule);
assert(max(abs([objval{:}])) == 0.0);

% Function for gradient computation
objVals = objective(model, states0, schedule);
objScaling = sum([objVals{:}]);
objectiveGradient = @(p) evalObjectiveBattmo(p, objective, simSetup, params, 'objScaling', objScaling);

%% Optional debug: Compare gradients using adjoints and finite difference approximation
debug = true;
if debug
    pTmp = getScaledParameterVector(simSetup, params);

    [vad, gad] = evalObjectiveBattmo(pTmp, objective, simSetup, params, ...
                                     'Gradient', 'AdjointAD');

    [vnum, gnum] = evalObjectiveBattmo(pTmp, objective, simSetup, params, ...
                                       'Gradient', 'PerturbationADNUM', 'PerturbationSize', 1e-7);
    [gad, gnum, abs(gad-gnum)]
end

%% Optimize
p0scaled = getScaledParameterVector(simSetup, params);
gradTol = 1e-7;
objChangeTol = 1e-4;
[v, pOptTmp, history] = unitBoxBFGS(p0scaled, objectiveGradient, ...
                                    'maximize', false, ...
                                    'gradTol', gradTol, ...
                                    'objChangeTol', objChangeTol, ...
                                    'logplot', true);
numIt = numel(history.val);

%% Unscale
pOpt = zeros(size(pOptTmp));
for k = 1:numel(pOpt)
    pOpt(k) = params{k}.unscale(pOptTmp(k));
end

% Compare with experimental params to see if we match
relErr = abs(pOpt - pExp) ./ pExp;
fprintf('pOpt             pExp            Rel err\n');
fprintf('%1.3g     \t %1.3g      \t %1.3g\n', [pOpt, pExp, relErr]');

%% Run model with pOpt params
jsonOpt = json;
for ip = 1:numel(params)
    loc = params{ip}.location;
    jsonOpt = setfield(jsonOpt, loc{:}, pOpt(ip));
end

[statesOpt, modelOpt] = runPImodel(jsonOpt);

%%
if do_plot
    set(0, 'defaultlinelinewidth', 2)

    getTime = @(states) cellfun(@(state) state.time, states);
    getE = @(states) cellfun(@(state) state.Control.E, states);

    t0 = getTime(states0);
    E0 = getE(states0);
    tOpt = getTime(statesOpt);
    EOpt = getE(statesOpt);
    tExp = getTime(statesExp);
    EExp = getE(statesExp);

    h = figure; hold on; grid on; axis tight
    plot(t0/hour, E0, 'displayname', 'E_{0}')
    plot(tExp/hour, EExp, '--', 'displayname', 'E_{exp}');
    plot(tOpt/hour, EOpt, 'displayname', 'E_{opt}')
    legend;

end

%% Summarize
pOrig = cellfun(@(p) p.getParameter(simSetup), params)';

fprintf('Initial guess:\n');
fprintf('%g\n', pOrig);
fprintf('\nOptimized values:\n');
fprintf('%g\n', pOpt);
fprintf('\nExperimental values:\n');
fprintf('%g\n', pExp);
fprintf('\nRelative error between optimized and experimental values:\n')
fprintf('%g\n', relErr);
fprintf('\nIterations:\n')
fprintf('%g\n', numIt);

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
