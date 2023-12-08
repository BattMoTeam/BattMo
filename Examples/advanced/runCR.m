%% Coin cell Lithium-Ion Battery Model
% This example demonstrates how to setup a simulation of a lithium
% manganese dioxide (CR) coin cell battery

%%

if mrstPlatform('octave')
    error('This demo cannot be run from Octave since Octave does not yet support the use of tables');
end

clear
close all

mrstVerbose off

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa upr

%% Define some shorthand names for simplicity.
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
co      = 'Coating';
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

compDims = table('rownames', {'NegativeCurrentCollector', ...
                              'NegativeCoating', ...
                              'Separator', ...
                              'PositiveCoating', ...
                              'PositiveCurrentCollector'});
numComponents = numel(compDims.Row);

%% Thickness
compDims.thickness = zeros(numComponents, 1);

compDims{'PositiveCoating', 'thickness'} = 67*micro*meter;
compDims{'Separator'      , 'thickness'} = 20*micro*meter;
compDims{'NegativeCoating', 'thickness'} = 50*micro*meter;

currentcollectors = {'PositiveCurrentCollector', 'NegativeCurrentCollector'};
compDims{currentcollectors, 'thickness'} = 0.5*(CRthickness - sum(compDims.thickness));

%% Diameters
compDims.diameter = zeros(numComponents, 1);

compDims{'PositiveCurrentCollector', 'diameter'} = 1;
compDims{'PositiveCoating'         , 'diameter'} = 0.8;
compDims{'Separator'               , 'diameter'} = 0.9;
compDims{'NegativeCoating'         , 'diameter'} = 0.8;
compDims{'NegativeCurrentCollector', 'diameter'} = 1;
compDims.diameter = compDims.diameter * CRdiameter;

%% Construct mesh
numRadial = 10;
numAngular = 50;
hz = min(compDims.thickness);
compDims.numCellLayers = max(2, round(compDims.thickness / hz));
compDims{currentcollectors, 'numCellLayers'} = ceil(round(0.25 * compDims{currentcollectors, 'numCellLayers'}));

disp(compDims);

params = struct('compDims'  , compDims , ...
                'numRadial' , numRadial, ...
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
Li_mass = masses.(ne).(co).val;

fprintf('Capacity %f mAh\n'  , C*1000/3600);
fprintf('Li content %f g\n'  , Li_mass * 1000);
fprintf('Battery mass %f g\n', mass * 1000);

%% Setup the time step schedule
% Smaller time steps are used to ramp up the current from zero to its
% operational value. Larger time steps are then used for the normal
% operation.
n         = 24 / 4;
dt        = [];
dt        = [dt; repmat(0.5e-4, n, 1).*1.5.^(1 : n)'];
totalTime = 2.0*hour/CRate;
n         = 20;
dt        = [dt; repmat(totalTime/n, n, 1)];
times     = [0; cumsum(dt)];
tt        = times(2 : end);
step      = struct('val', diff(times), 'control', ones(numel(tt), 1));

%% Setup the operating limits for the cell

control = model.Control.setupScheduleControl();

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
plotBatteryMesh(model, 'setstyle', false);

%% Run simulation and save output to folder
name = 'runCR';
dataFolder = 'BattMo';
problem = packSimulationProblem(initstate, model, schedule, dataFolder, 'Name', name, 'NonLinearSolver', nls);
problem.SimulatorSetup.OutputMinisteps = true;

clearSimulation = false;
if clearSimulation
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
figure
plot(time/hour, Inew, '-'); grid on; xlabel('time  / h');
figure
plot(time/hour, Enew, '-'); grid on; xlabel('time  / h');
if model.use_thermal
    % Plot half of the cell
    cells = model.(thermal).G.cells.centroids(:, 1) > 0;
    figure
    plotCellData(model.(thermal).G, states{end}.(thermal).T, cells);
    colorbar
    view(3), axis square tight
end

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
