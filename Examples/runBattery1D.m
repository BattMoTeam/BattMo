%% Battery 1D model
% 

% load MRST modules
mrstModule add ad-core mrst-gui mpfa

%%
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
thermal = 'ThermalModel';

%% Setup the computational mesh
% Here, we setup the 1D computational mesh that will be used for the
% simulation. The required discretization parameters are already included
% in the class BatteryGenerator1D. 
gen = BatteryGenerator1D();
gen.fac = 10;
gen = gen.applyResolutionFactors();

% Now, we update the paramobj with the properties of the mesh. 
paramobj = gen.updateBatteryInputParams(paramobj);

% !!! REMOVE THIS. SET THE RIGHT VALUES IN THE JSON !!! In this case, we
% change some of the values of the paramaters that were given in the json
% file to other values. This is done directly on the object paramobj. 
paramobj.(ne).(cc).EffectiveElectricalConductivity = 1e5;
paramobj.(pe).(cc).EffectiveElectricalConductivity = 1e5;
paramobj.(thermal).externalTemperature = paramobj.initT;

%%  Initialize the battery model. 
% The battery model is initialized by sending paramobj to the Battery class
% constructor. see :class:`Battery <Battery.Battery>`.
model = Battery(paramobj,'use_thermal',true,'use_solid_diffusion',true);
model.AutoDiffBackend= AutoDiffBackend();

%% Compute the nominal cell capacity and choose a discharge rate
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
fac=2;
total = 1.4*hour/CRate;
n=10;
dt0=total*1e-6;
times = getTimeSteps(dt0,n, total,fac);
dt= diff(times);
step = struct('val',diff(times),'control',ones(size(dt)));

% A stopping function is used to set the lower voltage cutoff limit.
pe = 'PositiveElectrode';
cc = 'CurrentCollector';
stopFunc = @(model, state, state_prev) (state.(pe).(cc).E < paramobj.Ucut); 

% A cource function is used to set the upper voltage cutoff limit. !!!
% Change this to an entry in the JSON with better variable names !!!
tup = 0.1; % rampup value for the current function, see rampupSwitchControl
inputE = 3.0; % Value when current control switches to voltage control
srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, inputI, inputE);

% we setup the control by assigning a source and stop function.
control = repmat(struct('src', srcfunc, 'stopFunction', stopFunc), 1, 1); 

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step); 

%% Setup the initial state
initstate = model.setupInitialState(); 

% Setup nonlinear solver 
nls = NonLinearSolver(); 
% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10; 
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false; 
nls.timeStepSelector=StateChangeTimeStepSelector('TargetProps',{{'PositiveElectrode','CurrentCollector','E'}},'targetChangeAbs',0.03);
% Change default tolerance for nonlinear solver
model.nonlinearTolerance = 1e-3*inputI; 
% Set verbosity
model.verbose = true;

% Run simulation
[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

%% Process output and recover the output voltage and current from the output states.
ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
Enew = cellfun(@(x) x.(pe).(cc).E, states); 
Inew = cellfun(@(x) x.(pe).(cc).I, states);
Tmax = cellfun(@(x) max(x.ThermalModel.T), states);
[SOCN,SOCP] =  cellfun(@(x) model.calculateSOC(x), states);
time = cellfun(@(x) x.time, states); 

%% Plot the the output voltage and current
plotDashboard(model,states,'step', 0)

% figure
% plot((time/hour), Enew, '*-', 'linewidth', 3)
% ylabel('Cell Voltage  /  V')
% xlabel('Time  /  h')
% setFigureStyle()
% 
% figure
% plot((time/hour), Inew, '*-', 'linewidth', 3)
% ylabel('Cell Current  /  A')
% xlabel('Time  /  h')
% setFigureStyle()
% 
% figure
% plot((time/hour), Tmax, '*-', 'linewidth', 3)
% ylabel('Max. Temperature  /  K')
% xlabel('Time  /  h')
% setFigureStyle()
% 
% figure
% plot((time/hour), [SOCP,SOCN], '*-', 'linewidth', 3)
% ylabel('SOC')
% xlabel('Time  /  h')
% legend('SOC positive','SOC negative')
% setFigureStyle()

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
