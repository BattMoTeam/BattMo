%% Setup
mrstModule add ad-core optimization mpfa

clear all
close all

do_plot = true;

%% Choose battery type

%json = parseBattmoJson(fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', 'lithium_ion_battery_nmc_graphite.json'));

json = parseBattmoJson(fullfile('ParameterData', 'ParameterSets', 'Xu2015', 'lfp.json'));

%% Generate "experimental" data
jsonExp = json;
statesExp = runPImodel(jsonExp);

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

%% Setup optimization
simulationSetup = struct('model', model, 'schedule', schedule, 'state0', state0);

% Convert config to ModelParameter (addParameter) input
parameters = configToModelParameter(simulationSetup, config);

% Objective function
objective = @(model, states, schedule, varargin) leastSquaresEI(model, states, schedule, 'statesRef', statesExp, varargin{:});

% Debug: the objective function evaluated at the experimental values should be zero
objval = objective(model, statesExp, schedule);
assert(max(abs([objval{:}])) == 0.0);

% Function for gradient computation
objVals = objective(model, states, schedule);
objScaling = sum([objVals{:}]);
objectiveGradient = @(p) evalMatchBattmo(p, objective, simulationSetup, parameters, 'objScaling', objScaling);

%% Debug: Compare gradients using adjoints and finite difference approximation
pTmp = getScaledParameterVector(simulationSetup, parameters);

[vad, gad] = evalMatchBattmo(pTmp, objective, simulationSetup, parameters, ...
                             'Gradient', 'AdjointAD');

[vnum, gnum] = evalMatchBattmo(pTmp, objective, simulationSetup, parameters, ...
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
relErr = abs(pOpt - pExp) ./ pExp
[pOpt, pExp, relErr]

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
fprintf('Initial guess\n');
p0'
fprintf('Optimized values\n');
pOpt
fprintf('Experimental values\n');
pExp
fprintf('Relative error between optimized and experimental values\n')
relErr
fprintf('Iterations\n')
numIt
