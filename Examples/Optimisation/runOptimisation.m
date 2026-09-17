clear all
close all

ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
am      = 'ActiveMaterial';
elyte   = 'Electrolyte';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
sep     = 'Separator';

jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_lco_graphite.json'));
jsonstruct.include_current_collectors = false;
jsonstruct.use_thermal = false;

jsonfilename = fullfile('Examples', 'JsonDataFiles', 'geometry1d.json');
jsonstruct_geometry = parseBattmoJson(jsonfilename);
jsonstruct = mergeStructs({jsonstruct, jsonstruct_geometry});

output = runBattery(jsonstruct, 'runSimulation', false, 'includeGridGenerator', true);

model   = output.model;
gridgen = output.gridGenerator;

css = CellSpecificationSummary(model);
css = css.setupPacking('packing_mass', 60*gram);

%% Setup default schedule

model.(ctrl).useCVswitch = true;
schedule = output.schedule;
model.(ctrl).rampupTime = 0.1;
schedule.control = model.(ctrl).setupScheduleControl();


%% Setup the initial state of the model

state0 = model.setupInitialState();

%% Setup the properties of the nonlinear solver

nls = NonLinearSolver();
nls.maxIterations  = 10;
nls.errorOnFailure = false;

model.nonlinearTolerance = 1e-5*model.Control.Imax;
model.verbose            = false;

dorunstartup = false;

if dorunstartup

    [wellSols, states, report] = simulateScheduleAD(state0, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); %#ok<UNRCH>

    ind = cellfun(@(x) not(isempty(x)), states);
    states = states(ind);
    E = cellfun(@(x) x.(ctrl).E, states);
    I = cellfun(@(x) x.(ctrl).I, states);
    time = cellfun(@(x) x.time, states);

    figure
    plot(time, E)

    figure
    plot(time, I)

end

%%

simulatorSetup = SimulationSetup(struct('model', model, ...
    'schedule', schedule, ...
    'initstate', state0, ...
    'NonLinearSolver', nls, ...
    'OutputMinisteps', false));

NPratio = 1.1;
paramsetter = NPlengthPorositySetter1D(model, gridgen, NPratio);

getValues = @(model, notused) paramsetter.getValues(model);
setValues = @(model, notused, v) paramsetter.setValues(model, v);

simulatorSetup.model = setValues(simulatorSetup.model, [], [80e-6; 0.3; 0.3]);
boxLims = [[50, 110]*1e-6; ...
           [0.15, 0.4]  ; ...
           [0.15, 0.4]];

% Rebuild the summary because its model is immutable, then restore the packing mass.
packingMass = css.packing_mass;
css = CellSpecificationSummary(simulatorSetup.model);
css = css.setupPacking('packing_mass', packingMass);
css.printSpecifications();

parameters{1} = ModelParameter(simulatorSetup            , ...
                               'name'     , 'lengthPoros', ...
                               'belongsTo', 'model'      , ...
                               'location' , {''}         , ...
                               'boxLims'  , boxLims      , ...
                               'getfun'   , getValues    , ...
                               'setfun'   , setValues);

objmatch = @(simsetup, states, varargin) EnergyOutput(simsetup, states, varargin{:});

options = {};

% setup objective function scaling. We use initial value
p = getScaledParameterVector(simulatorSetup, parameters);
v0 = evalObjectiveBattmo(p, objmatch, simulatorSetup, parameters, 'gradientMethod', 'None', options{:});

options = horzcat(options, {'objScaling', v0});


%%
doCompareGradient = false;

if doCompareGradient

    p = getScaledParameterVector(simulatorSetup, parameters); %#ok<UNRCH>
    [vad, gad]   = evalObjectiveBattmo(p, objmatch, simulatorSetup, parameters, 'gradientMethod', 'AdjointAD', options{:});
    perturbationSize = 1e-10;
    [vnum, gnum] = evalObjectiveBattmo(p, objmatch, simulatorSetup, parameters, ...
                                       'gradientMethod', 'PerturbationADNUM'  , ...
                                       'PerturbationSize', perturbationSize   , ...
                                       options{:});

    fprintf('Gradient computed using adjoint:\n');
    display(gad);
    fprintf('Numerical gradient:\n');
    display(gnum);

end

%%
doOptimization = true;

if doOptimization

    pBase = getScaledParameterVector(simulatorSetup, parameters);
    obj = @(p) evalObjectiveBattmo(p, objmatch, simulatorSetup, parameters, 'gradientMethod', 'AdjointAD', options{:});
    [v, pOpt, history] = unitBoxBFGS(pBase, obj, 'gradTol', 1e-4, 'objChangeTol', 1e-6);

    optimSetup = updateSetupFromScaledParameters(simulatorSetup, parameters, pOpt);

    doplot = true;

    if doplot

        [~, states] = simulateScheduleAD(simulatorSetup.initstate, simulatorSetup.model, simulatorSetup.schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);
        ind = cellfun(@(x) not(isempty(x)), states);
        states = states(ind);
        E = cellfun(@(x) x.(ctrl).E, states);
        time = cellfun(@(x) x.time, states);

        [~, states] = simulateScheduleAD(optimSetup.initstate, optimSetup.model, optimSetup.schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);
        ind = cellfun(@(x) not(isempty(x)), states);
        states = states(ind);
        EOpt = cellfun(@(x) x.(ctrl).E, states);
        timeOpt = cellfun(@(x) x.time, states);

        figure
        plot(time/hour, E, timeOpt/hour, EOpt)
        legend({'before', 'after'})
        grid on
        xlabel 'time  / h'
        ylabel 'E  / V'

    end

    vals    = parameters{1}.getParameter(simulatorSetup);
    valsOpt = parameters{1}.getParameter(optimSetup);

    % Print parameters
    fprintf('Initial and optimal parameter values (* indicates we hit boxLim)\n');
    for k = 1:numel(vals)
        name = paramsetter.fdnames{k};
        star = '';
        b = parameters{1}.boxLims;
        if abs(valsOpt(k)-b(k,1)) / b(k,1) < 1e-3 || abs(valsOpt(k)-b(k,2)) / b(k,2) < 1e-3
            star = '*';
        end

        fprintf('%s  \t %1.2e\t%1.2e %s\n', name, vals(k), valsOpt(k), star);
    end

    % Print statistics
    css.printSpecifications();
    cssOpt = CellSpecificationSummary(optimSetup.model);
    cssOpt = cssOpt.setupPacking('packing_mass', css.packing_mass);
    cssOpt.printSpecifications();

end



%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
