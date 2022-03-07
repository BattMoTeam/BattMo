%% Example: Pseudo-Three-Dimensional (P3D) Lithium-Ion Battery Model
% This example demonstrates how to setup a P3D model of a Li-ion battery
% and run a simple simulation.

% clear the workspace and close open figures
clear
close all
clc

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa

%% Setup the properties of Li-ion battery materials and cell design
% The properties and parameters of the battery cell, including the
% architecture and materials, are set using an instance of
% :class:`BatteryInputParams <Battery.BatteryInputParams>`. This class is
% used to initialize the simulation and it propagates all the parameters
% throughout the submodels. The input parameters can be set manually or
% provided in json format. All the parameters for the model are stored in
% the paramobj object.
jsonstruct = parseBatmoJson('JsonDatas/lithiumbattery.json');
paramobj = BatteryInputParams(jsonstruct);

% We define some shorthand names for simplicity.
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
eac     = 'ElectrodeActiveComponent';
cc      = 'CurrentCollector';
elyte   = 'Electrolyte';
sep     = 'Separator';
thermal = 'ThermalModel';

%% Setup the geometry and computational mesh
% Here, we setup the 2D computational mesh that will be used for the
% simulation. The required discretization parameters are already included
% in the class BatteryGenerator2D. 
gen = BatteryGenerator2D();

% Now, we update the paramobj with the properties of the mesh.
paramobj = gen.updateBatteryInputParams(paramobj);

% !!! REMOVE THIS. SET THE RIGHT VALUES IN THE JSON !!! In this case, we
% change some of the values of the paramaters that were given in the json
% file to other values. This is done directly on the object paramobj. 
paramobj.(ne).(cc).EffectiveElectricalConductivity = 1e5;
paramobj.(pe).(cc).EffectiveElectricalConductivity = 1e5;

%%  Initialize the battery model. 
% The battery model is initialized by sending paramobj to the Battery class
% constructor. see :class:`Battery <Battery.Battery>`.
model = Battery(paramobj);

%% Plot the mesh
% The mesh is plotted using the plotGrid() function from MRST. 
figure
plotGrid(model.(ne).(cc).G,     'facecolor', [191, 124, 84]./255,   'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
plotGrid(model.(ne).(eac).G,    'facecolor', [13, 13, 13]./255,     'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
plotGrid(model.(elyte).(sep).G, 'facecolor', [239, 204, 97]./255,   'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
plotGrid(model.(pe).(eac).G,    'facecolor', [64, 64, 64]./255,     'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
plotGrid(model.(pe).(cc).G,     'facecolor', [191, 191, 191]./255,  'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
axis tight;
legend({'negative electrode current collector' , ...
        'negative electrode active material'   , ...
        'separator'                            , ...
        'positive electrode active material'   , ...
        'positive electrode current collector'}, ...
       'location', 'south west'),
setFigureStyle('quantity', 'single');
drawnow();
pause(0.1);

%% Compute the nominal cell capacity and choose a C-Rate
% The nominal capacity of the cell is calculated from the active materials.
% This value is then combined with the user-defined C-Rate to set the cell
% operational current. 
C      = computeCellCapacity(model);
CRate  = 1;
inputI = (C/hour)*CRate; % current 

%% Setup the time step schedule 
% Smaller time steps are used to ramp up the current from zero to its
% operational value. Larger time steps are then used for the normal
% operation. 
n           = 25; 
dt          = []; 
dt          = [dt; repmat(0.5e-4, n, 1).*1.5.^[1:n]']; 
totalTime   = 1.4*hour/CRate;
n           = 40; 
dt          = [dt; repmat(totalTime/n, n, 1)]; 
times       = [0; cumsum(dt)]; 
tt          = times(2 : end); 
step        = struct('val', diff(times), 'control', ones(numel(tt), 1)); 

%% Setup the operating limits for the cell
% The maximum and minimum voltage limits for the cell are defined using
% stopping and source functions. A stopping function is used to set the
% lower voltage cutoff limit. A source function is used to set the upper
% voltage cutoff limit. 
stopFunc    = @(model, state, state_prev) (state.(pe).(cc).E < 2.0); 
tup         = 0.1; % rampup value for the current function, see rampupSwitchControl
inputE      = 3; % Value when current control switches to voltage control
srcfunc     = @(time, I, E) rampupSwitchControl(time, tup, I, E, inputI, inputE);

% we setup the control by assigning a source and stop function.
control = repmat(struct('src', srcfunc, 'stopFunction', stopFunc), 1, 1); 

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step); 

%% Setup the initial state of the model
% The initial state of the model is dispatched using the
% model.setupInitialState()method. 
initstate = model.setupInitialState(); 

%% Setup the properties of the nonlinear solver 
nls = NonLinearSolver(); 
% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10; 
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false; 
% Timestep selector
nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', ...
                                                  {{'PositiveElectrode', 'CurrentCollector', 'E'}}, ...
                                                  'targetChangeAbs', 0.03);

% Change default tolerance for nonlinear solver
model.nonlinearTolerance = 1e-5; 
% Set verbosity of the solver (if true, value of the residuals for every equation is given)
model.verbose = true;

%% Run simulation
[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, ...
                                                'OutputMinisteps', true, ...
                                                'NonLinearSolver', nls); 

%%  Process output and recover the output voltage and current from the output states.
ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
Enew = cellfun(@(x) x.(pe).(cc).E, states); 
Inew = cellfun(@(x) x.(pe).(cc).I, states);
time = cellfun(@(x) x.time, states); 

%% Plot an animated summary of the results
plotDashboard(model, states, 'step', 0);

%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BatMo

BatMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BatMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BatMo.  If not, see <http://www.gnu.org/licenses/>.
%}
