if mrstPlatform('octave')
    error('This demo cannot be run from Octave since Octave does not yet support the use of tables');
end

%% Setup
mrstModule add ad-core optimization mpfa

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

%% Choose battery type

jsonParams = parseBattmoJson(fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', 'lithium_ion_battery_nmc_graphite.json'));

% jsonParams = parseBattmoJson(fullfile('ParameterData', 'ParameterSets', 'Xu2015', 'lfp.json'));
% jsonGeom = parseBattmoJson(fullfile('Examples', 'jsondatafiles', 'geometry1d.json'));
% json = mergeJsonStructs({jsonParams, jsonGeom});

json = jsonParams;

%% Generate "experimental" data
jsonExp = json;
statesExp = runPImodel(jsonExp);
%outputExp = runBatteryJson(jsonExp);

%% Setup parameter identification problem

% Get configuration for matching Butler-Volmer k0 parameters
%config = configButlerVolmer(jsonExp);

% Get configuration for matching Bruggeman coefficients
%config = configBruggeman(jsonExp);

% Get configuration for matching volumetric surface areas
%config = configVolumetricSurfaceArea(jsonExp);

% Get configuration for parameter identification of several parameters
% (Butler-Volmer k0 parameters, Bruggeman coefficients and volumetric
% surface areas)
config = configMultiple(jsonExp);

%% Run with initial guess
json0 = json;
numVars = numel(config.Row);
for k = 1:numVars
    loc = config.location{k};
    json0 = setfield(json0, loc{:}, config.initialGuess(k));
end
[states, model, schedule, state0] = runPImodel(json0);
% output = runBatteryJson(json0);

%% Setup optimization
simulationSetup = struct('model', model, 'schedule', schedule, 'state0', state0);

% Convert config to ModelParameter (addParameter) input
parameters = configToModelParameter(simulationSetup, config);

% Objective function
objective = @(model, states, schedule, varargin) leastSquaresEI(model, states, statesExp, schedule, varargin{:});

% Debug: the objective function evaluated at the experimental values should be zero
objval = objective(model, statesExp, schedule);
assert(max(abs([objval{:}])) == 0.0);

% Function for gradient computation
objVals = objective(model, states, schedule);
objScaling = sum([objVals{:}]);
objectiveGradient = @(p) evalObjectiveBattmo(p, objective, simulationSetup, parameters, 'objScaling', objScaling);

%% Debug: Compare gradients using adjoints and finite difference approximation
pTmp = getScaledParameterVector(simulationSetup, parameters);

[vad, gad] = evalObjectiveBattmo(pTmp, objective, simulationSetup, parameters, ...
                             'Gradient', 'AdjointAD');

[vnum, gnum] = evalObjectiveBattmo(pTmp, objective, simulationSetup, parameters, ...
                               'Gradient', 'PerturbationADNUM', 'PerturbationSize', 1e-7);
[gad, gnum]

%% Optimize
p0scaled = getScaledParameterVector(simulationSetup, parameters);
gradTol = 1e-7;
objChangeTol = 1e-4;
[v, pOptTmp, history] = unitBoxBFGS(p0scaled, objectiveGradient, 'maximize', false, ...
                                    'gradTol', gradTol, ...
                                    'objChangeTol', objChangeTol);
numIt = numel(history.val);

%% Unscale
pOpt = zeros(size(pOptTmp));
for k = 1:numel(pOpt)
    pOpt(k) = parameters{k}.unscale(pOptTmp(k));
end

% Compare with reference parameters
pExp = cellfun(@(x) x.referenceValue, parameters);
relErr = abs(pOpt - pExp) ./ pExp;
fprintf('%g \t %g \t %g\n', [pOpt, pExp, relErr]');

%% Run model with pOpt parameters
jsonOpt = json;
for k = 1:numVars
    loc = parameters{k}.location;
    jsonOpt = setfield(jsonOpt, loc{:}, pOpt(k));
end

[statesOpt, modelOpt] = runPImodel(jsonOpt);
p0 = cellfun(@(x) x.getParameterValue(simulationSetup), parameters)';

%%
if do_plot

    [t0, E0] = extractTE(states);
    [tOpt, EOpt] = extractTE(statesOpt);
    [tExp, EExp] = extractTE(statesExp);

    h = figure;
    plot(t0, E0, tExp, EExp, '--', tOpt, EOpt, 'linewidth', 2)
    grid on

    legstr = @(x, s) [s, sprintf(repmat(' %1.2g', 1, numel(x)), x)];
    legend(legstr(p0, 'E_{0}'), ...
           legstr(pExp, 'E_{ref}'), ...
           legstr(pOpt, 'E_{opt}'), ...
           'location', 'SW');

    tit = sprintf(['relErr=[', repmat('%1.1g ', 1, numVars), '] it=%d'], ...
                  relErr, numel(history.u));
    title(tit)

end

%% Summarize
fprintf('Initial guess:\n');
fprintf('%g\n', p0);
fprintf('Optimized values:\n');
fprintf('%g\n', pOpt);
fprintf('Experimental values:\n');
fprintf('%g\n', pExp);
fprintf('Relative error between optimized and experimental values:\n')
fprintf('%g\n', relErr);
fprintf('Iterations:\n')
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
