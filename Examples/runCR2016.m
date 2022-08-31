clear all
close all
clc

mrstVerbose off

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa upr

mrstDebug(99)

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
jsonstruct = parseBattmoJson('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json');

jsonstruct.use_thermal = false;

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
CR2016_diameter = 20*milli*meter;
CR2016_thickness = 1.6*milli*meter;

compnames = {'NegativeCurrentCollector', ...
             'NegativeActiveMaterial', ...
             'ElectrolyteSeparator', ...
             'PositiveActiveMaterial', ...
             'PositiveCurrentCollector'};
num_components = 5;
compdims = table('rownames', compnames);

%% Thickness
% From Joule paper: this gives only half of the Li content compared to
% Energizer spread sheet, and 1/20 of the typical capacity. Hence the factor 2.
thickness_factor = 2;
compdims.thickness = zeros(num_components, 1);

compdims{'PositiveActiveMaterial', 'thickness'} = 67*micro*meter * thickness_factor;
compdims{'ElectrolyteSeparator'  , 'thickness'} = 20*micro*meter * thickness_factor;
compdims{'NegativeActiveMaterial', 'thickness'} = 50*micro*meter * thickness_factor;

ccs = {'PositiveCurrentCollector', 'NegativeCurrentCollector'};
compdims{ccs, 'thickness'} = 0.5*(CR2016_thickness - sum(compdims.thickness));

%% Diameters
compdims.diameter = zeros(num_components, 1);

compdims{'PositiveCurrentCollector', 'diameter'} = 1;
compdims{'PositiveActiveMaterial'  , 'diameter'} = 0.9;
compdims{'ElectrolyteSeparator'    , 'diameter'} = 0.7;
compdims{'NegativeActiveMaterial'  , 'diameter'} = 0.9;
compdims{'NegativeCurrentCollector', 'diameter'} = 1;
compdims.diameter = compdims.diameter * CR2016_diameter;

%% 
meshSize = max(compdims.diameter) / 10;
hz = min(compdims.thickness);
compdims.numCellLayers = max(2, round(compdims.thickness / hz));
compdims{ccs, 'numCellLayers'} = ceil(round(0.25 * compdims{ccs, 'numCellLayers'}));

%use_sector = true;
use_sector = false;

offset = [0.0, 0.0];

params = struct('compdims', compdims, ...
                'meshSize', meshSize, ...
                'use_sector', use_sector, ...
                'offset', offset);

gen = CoinCellBatteryGenerator();

% Now, we update the paramobj with the properties of the mesh.
paramobj = gen.updateBatteryInputParams(paramobj, params);

%%  Initialize the battery model.
% The battery model is initialized by sending paramobj to the Battery class
% constructor. see :class:`Battery <Battery.Battery>`.
model = Battery(paramobj);

% % test
% idx = isnan(model.(elyte).volumeFraction);
% model.(elyte).volumeFraction(idx) = 0.5;


%% Compute the nominal cell capacity and choose a C-Rate
% The nominal capacity of the cell is calculated from the active materials.
% This value is then combined with the user-defined C-Rate to set the cell
% operational current.
C     = computeCellCapacity(model);
CRate = 1;

fprintf('Capacity %f mAh\n', C*1000/3600)

[mass, masses] = computeCellMass(model);
fprintf('Li content %f g\n', masses.(ne).(am).val*1000);
fprintf('Battery mass %f g\n', mass*1000);



%% Setup the time step schedule
% Smaller time steps are used to ramp up the current from zero to its
% operational value. Larger time steps are then used for the normal
% operation.
n         = 24 / 4;
dt        = [];
dt        = [dt; repmat(0.5e-4, n, 1).*1.5.^[1 : n]'];
totalTime = 2.0*hour/CRate;
n         = 20; %40;
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
plotGrid(model.(elyte).G,model.(elyte).G.cells.centroids(:, 1)>0,       'facecolor', colors(6,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1], 'facealpha', 0.1);view(3);title('elyte')

drawnow

return



% %% Run simulation
% mrstVerbose off
% [wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

if ~exist('tag')
    tag = datetime;
end
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
plot(time, Inew, '-'), title('I', num2str(diamfactor)); grid on; xlabel 'time (s)'
figure
plot(time, Enew, '-'), title('E', num2str(diamfactor)); grid on; xlabel 'time (s)'

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
