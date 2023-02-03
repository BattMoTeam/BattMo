clear all
close all
clc

mrstVerbose off

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa upr

%% Define some shorthand names for simplicity.
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
eac     = 'ActiveMaterial';
cc      = 'CurrentCollector';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
ctrl    = 'Control';
am      = 'ActiveMaterial';
sep     = 'Separator';

%% Setup the properties of Li-ion battery materials and cell design
jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
jsonstruct.use_thermal = false;
jsonstruct.include_current_collectors = true;

paramobj = BatteryInputParams(jsonstruct);

use_cccv = false;
if use_cccv
    cccvstruct = struct('controlPolicy'     , 'CCCV',  ...
                        'CRate'             , 1         , ...
                        'lowerCutoffVoltage', 2         , ...
                        'upperCutoffVoltage', 4.1       , ...
                        'dIdtLimit'         , 0.01      , ...
                        'dEdtLimit'         , 0.01);
    cccvparamobj = CcCvControlModelInputParams(cccvstruct);
    paramobj.Control = cccvparamobj;
end

%% Setup the geometry and grid for the components
CRdiameter = 20*milli*meter;
CRthickness = 1.6*milli*meter;

compNames = {'NegativeCurrentCollector', ...
             'NegativeActiveMaterial', ...
             'ElectrolyteSeparator', ...
             'PositiveActiveMaterial', ...
             'PositiveCurrentCollector'};
numComponents = numel(compNames);
compDims = table('rownames', compNames);

%% Thickness
compDims.thickness = zeros(numComponents, 1);

compDims{'PositiveActiveMaterial', 'thickness'} = 67*micro*meter;
compDims{'ElectrolyteSeparator'  , 'thickness'} = 20*micro*meter;
compDims{'NegativeActiveMaterial', 'thickness'} = 50*micro*meter;

currentcollectors = {'PositiveCurrentCollector', 'NegativeCurrentCollector'};
compDims{currentcollectors, 'thickness'} = 0.5*(CRthickness - sum(compDims.thickness));

%% Diameters
compDims.diameter = zeros(numComponents, 1);

compDims{'PositiveCurrentCollector', 'diameter'} = 1;
compDims{'PositiveActiveMaterial'  , 'diameter'} = 0.8;
compDims{'ElectrolyteSeparator'    , 'diameter'} = 0.9;
compDims{'NegativeActiveMaterial'  , 'diameter'} = 0.8;
compDims{'NegativeCurrentCollector', 'diameter'} = 1;
compDims.diameter = compDims.diameter * CRdiameter;

%% Construct mesh
numRadial = 10;
numAngular = 50;
hz = min(compDims.thickness);
compDims.numCellLayers = max(2, round(compDims.thickness / hz));
compDims{currentcollectors, 'numCellLayers'} = ceil(round(0.25 * compDims{currentcollectors, 'numCellLayers'}));

params = struct('compDims', compDims, ...
                'numRadial', numRadial, ...
                'numAngular', numAngular);

gen = CoinCellBatteryGenerator();

% Now, we update the paramobj with the properties of the mesh.
paramobj = gen.updateBatteryInputParams(paramobj, params);

%%  Initialize the battery model.
% The battery model is initialized by sending paramobj to the Battery class
% constructor. see :class:`Battery <Battery.Battery>`.
model = Battery(paramobj);

%% Compute the nominal cell capacity and choose a C-Rate
% The nominal capacity of the cell is calculated from the active materials.
% This value is then combined with the user-defined C-Rate to set the cell
% operational current.
C = computeCellCapacity(model);
CRate = 1;

% Compute masses
[mass, masses] = computeCellMass(model);
Li_mass = masses.(ne).(am).val;

fprintf('Capacity %f mAh\n', C*1000/3600);
fprintf('Li content %f g\n', Li_mass * 1000);
fprintf('Battery mass %f g\n', mass * 1000);

%% Setup the time step schedule
% Smaller time steps are used to ramp up the current from zero to its
% operational value. Larger time steps are then used for the normal
% operation.
n         = 24 / 4;
dt        = [];
dt        = [dt; repmat(0.5e-4, n, 1).*1.5.^[1 : n]'];
totalTime = 2.0*hour/CRate;
n         = 20;
dt        = [dt; repmat(totalTime/n, n, 1)];
times     = [0; cumsum(dt)];
tt        = times(2 : end);
step      = struct('val', diff(times), 'control', ones(numel(tt), 1));

%% Setup the operating limits for the cell

switch model.Control.controlPolicy
  case 'IEswitch'
    % The maximum and minimum voltage limits for the cell are defined using
    % stopping and source functions. A stopping function is used to set the
    % lower voltage cutoff limit. A source function is used to set the upper
    % voltage cutoff limit.
    tup = 0.1; % rampup value for the current function, see rampupSwitchControl
    srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                model.Control.Imax, ...
                                                model.Control.lowerCutoffVoltage);
    % setup the control by assigning a source and stop function.
    control = struct('src', srcfunc, 'IEswitch', true);

  case 'CCCV'
    control = struct('CCCV', true);

  otherwise
    error('control policy not recognized');
end

    
% This control is used to set up the schedule
schedule = struct('control', control, 'step', step);

%% Setup the initial state of the model
% The initial state of the model is dispatched using the
% model.setupInitialState()method.
initstate = model.setupInitialState();

initstate.(ctrl).I = 0;
initstate.(ctrl).ctrlType = 'constantCurrent';

%% Setup the properties of the nonlinear solver
nls = NonLinearSolver();
nls.maxIterations = 10;
nls.errorOnFailure = false;
nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', ...
                                                   {{ctrl, 'E'}}, ...
                                                   'targetChangeAbs', 0.03);
model.nonlinearTolerance = 1e-5;
model.verbose = true;

%% Plot
colors = crameri('vik', 6);
figure
plotGrid(model.(ne).(cc).G,     'facecolor', colors(1,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
plotGrid(model.(ne).(am).G,     'facecolor', colors(2,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
plotGrid(model.(elyte).(sep).G, 'facecolor', colors(3,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
plotGrid(model.(pe).(am).G,     'facecolor', colors(4,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
plotGrid(model.(pe).(cc).G,     'facecolor', colors(5,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
axis tight;
legend({'negative electrode current collector' , ...
        'negative electrode active material'   , ...
        'separator'                            , ...
        'positive electrode active material'   , ...
        'positive electrode current collector'}, ...
       'location', 'southwest')
view(-37, 14)
drawnow


%% Additional plots for visualizing the different components
doplot = false;
if doplot
    figure
    plotGrid(model.(ne).(cc).G,     'facecolor', colors(1,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);view(3);title('ne cc')
    figure
    plotGrid(model.(ne).(am).G,     'facecolor', colors(2,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);view(3);title('ne am');
    figure
    plotGrid(model.(elyte).(sep).G, 'facecolor', colors(3,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);view(3);title('elyte sep')
    figure
    plotGrid(model.(pe).(am).G,     'facecolor', colors(4,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);view(3);title('pe am')
    figure
    plotGrid(model.(pe).(cc).G,     'facecolor', colors(5,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);view(3);title('pe cc')
    figure
    plotGrid(model.(elyte).G,model.(elyte).G.cells.centroids(:, 1)>0,       'facecolor', colors(6,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1], 'facealpha', 0.1);view(3);title('elyte')
    drawnow
end

%% Run simulation and save output to folder
if ~exist('tag')
    tag = datetime;
end
name = sprintf('cr2016-%s', tag);
dataFolder = sprintf('BattMo_%s', date);
problem = packSimulationProblem(initstate, model, schedule, dataFolder, 'Name', name, 'NonLinearSolver', nls);
problem.SimulatorSetup.OutputMinisteps = true;

resetSimulation = false
if resetSimulation
    % Clear previously computed simulation
    clearPackedSimulatorOutput(problem, 'prompt', false);
end
simulatePackedProblem(problem);
[globvars, states, report] = getPackedSimulatorOutput(problem);


%%  Process output and recover the output voltage and current from the output states.
ind = cellfun(@(x) not(isempty(x)), states);
states = states(ind);
Enew = cellfun(@(x) x.(ctrl).E, states);
Inew = cellfun(@(x) x.(ctrl).I, states);
time = cellfun(@(x) x.time, states);

%% Plot results
figure, hold on
plot(time, Inew, '-'); title('I'); grid on; xlabel 'time (s)'
figure, hold on
plot(time, Enew, '-'); title('E'); grid on; xlabel 'time (s)'
if model.use_thermal
    figure
    plotCellData(model.(thermal).G, states{end}.(thermal).T, model.(thermal).G.cells.centroids(:, 1) > 0)
    colorbar
end

%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
