clear all
close all
clc

mrstVerbose off

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa upr

%% Setup the properties of Li-ion battery materials and cell design
jsonstruct = parseBattmoJson('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json');
paramobj = BatteryInputParams(jsonstruct);

% We define some shorthand names for simplicity.
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
eac     = 'ActiveMaterial';
cc      = 'CurrentCollector';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
ctrl    = 'Control';
am      = 'ActiveMaterial';
sep     = 'Separator';

%% Setup the geometry and computational mesh
CR2016_diameter = 20*milli*meter;
CR2016_thickness = 1.6*milli*meter;
%CR2016_thickness = 0.25*milli*meter;

components = {'NegativeCurrentCollector', ...
              'NegativeActiveMaterial', ...
              'ElectrolyteSeparator', ...
              'PositiveActiveMaterial', ...
              'PositiveCurrentCollector'};

% zlength = [10; 100; 50; 100; 10];
% zz = CR2016_thickness * zlength / sum(zlength);
% thickness = containers.Map(components, zz);

thickness = containers.Map();
thickness('PositiveCurrentCollector') = 0.1*milli*meter;
thickness('PositiveActiveMaterial')   = 0.45*milli*meter;
thickness('ElectrolyteSeparator')     = 0.25*milli*meter;
thickness('NegativeActiveMaterial')   = 0.7*milli*meter;
thickness('NegativeCurrentCollector') = thickness('PositiveCurrentCollector');

tag = '0.6-0.7-realz'
ddfactor = [1, 0.6, 0.7, 0.6, 1];
% tag = 'sep0.9-am0.8';
% ddfactor = [1, 0.8, 0.9, 0.8, 1];
dd = ddfactor * CR2016_diameter;
%dd = ones(1, 5) * CR2016_diameter;
diameter = containers.Map(components, dd);


meshSize = max(cell2mat(diameter.values)) / 10;

numCellLayers = containers.Map();
hz = 3*min(cell2mat(thickness.values));
for k = keys(thickness)
    key = k{1};
    numCellLayers(key) = max(2, round(thickness(key)/hz));
end

use_sector = true;
%use_sector = false;

params = struct('thickness', thickness, ...
                'diameter', diameter, ...
                'numCellLayers', numCellLayers, ...
                'meshSize', meshSize, ...
                'use_sector', use_sector);

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
C     = computeCellCapacity(model);
CRate = 1;

fprintf('Capacity %f mAh\n', C*1000/3600)

[~, masses] = computeCellMass(model);
fprintf('Li content %f g\n', masses.(ne).(am).val*1000);



%% Setup the time step schedule
% Smaller time steps are used to ramp up the current from zero to its
% operational value. Larger time steps are then used for the normal
% operation.
n         = 24 / 4;
dt        = [];
dt        = [dt; repmat(0.5e-4, n, 1).*1.5.^[1 : n]'];
totalTime = 60*second /CRate; %1.4*hour/CRate;
n         = 20 %40;
dt        = [dt; repmat(totalTime/n, n, 1)];
times     = [0; cumsum(dt)];
tt        = times(2 : end);
step      = struct('val', diff(times), 'control', ones(numel(tt), 1));

%% Setup the operating limits for the cell
% The maximum and minimum voltage limits for the cell are defined using
% stopping and source functions. A stopping function is used to set the
% lower voltage cutoff limit. A source function is used to set the upper
% voltage cutoff limit.
tup = 0.1; % rampup value for the current function, see rampupSwitchControl
srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                            model.Control.Imax, ...
                                            model.Control.lowerCutoffVoltage);
% we setup the control by assigning a source and stop function.
control = struct('src', srcfunc, 'IEswitch', true);

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step);

%% Setup the initial state of the model
% The initial state of the model is dispatched using the
% model.setupInitialState()method.
initstate = model.setupInitialState();

%% Setup the properties of the nonlinear solver
nls = NonLinearSolver();
nls.maxIterations = 10;
nls.errorOnFailure = false;
nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', ...
                                                  {{ctrl, 'E'}}, ...
                                                   'targetChangeAbs', 0.03);

mrstModule add linearsolvers
%nls.LinearSolver = AMGCLSolverAD('verbose', true, 'reduceToCell', false);
%nls.LinearSolver = AGMGSolverAD('verbose', true, 'reduceToCell', true);
% nls.LinearSolver.tolerance = 1e-3;
% nls.LinearSolver.maxIterations = 30;
nls.maxIterations = 10;
nls.verbose = 10;



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
%plotGrid(model.(elyte).G,       'facecolor', colors(6,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1], 'facealpha', 0.1);
axis tight;
legend({'negative electrode current collector' , ...
        'negative electrode active material'   , ...
        'separator'                            , ...
        'positive electrode active material'   , ...
        'positive electrode current collector'}, ...
       'location', 'southwest')
view(3)
drawnow

%return

%% More plot
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
plotGrid(model.(elyte).G,       'facecolor', colors(6,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1], 'facealpha', 0.1);view(3);title('elyte')

%return



% %% Run simulation
% mrstVerbose off
% [wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);


name = sprintf('cr2016-%s', tag);
dataFolder = sprintf('BattMo_%s', date);
problem = packSimulationProblem(initstate, model, schedule, dataFolder, 'Name', name, 'NonLinearSolver', nls);
problem.SimulatorSetup.OutputMinisteps = true;

resetSimulation = true;
if resetSimulation
    %% clear previously computed simulation
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

%% Plot an animated summary of the results
%plotDashboard(model, states, 'step', 0);


figure
plot(time, Inew, '.-'), title('I', ddfactor); grid on
figure
plot(time, Enew, '.-'), title('E', ddfactor); grid on

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
