if mrstPlatform('octave')
    if compare_versions(version, "8.1", "<=")
        error('This demo cannot be run in Octave versions less than 8.1, since the ind2sub behaves differently from MATLAB');
    end
end

% Setup mrst modules

mrstModule add ad-core mrst-gui mpfa agmg

%% We setup the geometrical parameters for a 4680 battery.
%% Those will be gathered in structure spiralparams (see below) and used by SpiralBatteryGenerator to generate the spiral layered geometry of the jelly roll

% widths of each component ordered as
% - positive current collector
% - positive electrode
% - electrolyte separator
% - negative electrode
% - negative current collector

widths = [25, 64, 15, 57, 15]*micro*meter;

widthDict = containers.Map(...
    {'Separator'               , ...
     'NegativeCoating'         , ...
     'NegativeCurrentCollector', ...
     'PositiveCoating'         , ...
     'PositiveCurrentCollector'}, ...
    widths);

nwidths = [widthDict('PositiveCoating')          ;...
           widthDict('PositiveCurrentCollector') ;...
           widthDict('PositiveCoating')          ;...
           widthDict('Separator')                ;...
           widthDict('NegativeCoating')          ;...
           widthDict('NegativeCurrentCollector') ;...
           widthDict('NegativeCoating')          ;...
           widthDict('Separator')];

dr = sum(nwidths);

% Radii of the jelly roll
rInner = 2*milli*meter;
rOuter = 46*milli*meter/2;
dR = rOuter - rInner;

% Height of the jelly roll
L = 80*milli*meter;

% Computed number of windings
nwindings = ceil(dR/dr);

% Number of discretization cells in radial direction for each component.
nrDict = containers.Map({'Separator'               , ...
                         'NegativeCoating'         , ...
                         'NegativeCurrentCollector', ...
                         'PositiveCoating'         , ...
                         'PositiveCurrentCollector'}, ...
                        [3, 3, 3, 3, 3]);

% Number of discretization cells in the longitudonal
nL  = 5;
nas = 3;

% Structure that describes the tab setups (see SpiralBatteryGenerator)
tabparams.tabcase   = 'aligned tabs';
tabparams.width     = 3*milli*meter;
tabparams.fractions = linspace(0.01, 0.9, 6);

spiralparams = struct('nwindings'   , nwindings, ...
                      'rInner'      , rInner   , ...
                      'widthDict'   , widthDict, ...
                      'nrDict'      , nrDict   , ...
                      'nas'         , nas      , ...
                      'L'           , L        , ...
                      'nL'          , nL       , ...
                      'tabparams'   , [], ...
                      'angleuniform', true);

% The input material parameters given in json format are used to populate the paramobj object.
jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
jsonstruct.include_current_collectors = true;

paramobj = BatteryInputParams(jsonstruct);

th = 'ThermalModel';
%paramobj.(th).externalHeatTransferCoefficientSideFaces = 100*watt/meter^2;
%paramobj.(th).externalHeatTransferCoefficientTopFaces = 10*watt/meter^2;
paramobj.(th).externalHeatTransferCoefficient = 10;

gen = SectorBatteryGenerator();

paramobj = gen.updateBatteryInputParams(paramobj, spiralparams);

model = Battery(paramobj);

%% Setup schedule
CRate = 1;
fac   = 2;
total = 1.4*hour/CRate;
n     = 10;
dt0   = total*1e-6;
times = getTimeSteps(dt0, n, total, fac);

% Compute the cell capacity, which used to compute schedule from CRate
capacity = computeCellCapacity(model);
inputI   = (capacity/hour)*CRate;
inputE   = 3;

tt   = times(2 : end);
step = struct('val', diff(times), 'control', ones(numel(tt), 1));
tup  = 0.1/CRate;

simcase = 'discharge';

switch simcase

  case 'discharge'
    srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                model.Control.Imax, ...
                                                model.Control.lowerCutoffVoltage);

    % Setup the control by assigning a source and stop function.
    control = struct('src', srcfunc, 'IEswitch', true);
    schedule  = struct('control', control, 'step', step);

    %% We setup the initial state
    initstate = model.setupInitialState();
    initElytec = 1*mol/litre;
    initstate.Electrolyte.c = initElytec*ones(model.Electrolyte.G.cells.num, 1);

  case 'charge'

    model.SOC = 0.01;
    initstate = model.setupInitialState();
    srcfunc  = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                 - model.Control.Imax, ...
                                                 model.Control.lowerCutoffVoltage);
    control = struct('src', srcfunc, 'IEswitch', true);
    schedule = struct('control', control, 'step', step);

  otherwise
    error('simcase not recognized')

end

% Setup nonlinear solver
nls = NonLinearSolver();

% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10;
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false;
% Change default tolerance for nonlinear solver
model.nonlinearTolerance = 1e-4;
% Get more or less verbose output
model.verbose = true;

% Options for experimenting with linear solver
nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', {{'Control', 'E'}}, 'targetChangeAbs', 0.03);

use_amg = false;
if use_amg
    mrstModule add agmg
    nls.LinearSolver = LinearSolverBattery('method', 'agmg', 'verbosity', 0);
    nls.maxIterations = 10;
end

% Run simulation
dataFolder = 'BattMo';
problem = packSimulationProblem(initstate, model, schedule, dataFolder, 'Name', 'jellyroll', 'NonLinearSolver', nls);
problem.SimulatorSetup.OutputMinisteps = true;

clearSimulation = true;
if clearSimulation
    % clear previously computed simulation
    clearPackedSimulatorOutput(problem, 'prompt', false);
end
simulatePackedProblem(problem);
[globvars, states, report] = getPackedSimulatorOutput(problem);


%% Plot states
figure
plotToolbar(model.G, states);
view([0,-1,0]);

%%  Process output and recover the output voltage and current from the output states.

ind = cellfun(@(x) not(isempty(x)), states);
states = states(ind);
ctrl = 'Control';
E = cellfun(@(x) x.(ctrl).E, states);
I = cellfun(@(x) x.(ctrl).I, states);
time = cellfun(@(x) x.time, states);

figure
plot(time, E);

%% Plot an animated summary of the results
%plotDashboard(model, states, 'step', 0);



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
