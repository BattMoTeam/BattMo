%% Battery 1D model
% 

% load MRST modules
mrstModule add ad-core mrst-gui mpfa

%%
% We create an instance of :class:`BatteryInputParams <Battery.BatteryInputParams>`. This class is used to initiate the
% battery simulator and it propagates all the parameters through out the submodels.

% The input parameters can be given in json format. The json file is read and used to populate the paramobj object.
jsonstruct = parseBatmoJson('JsonDatas/lithiumbattery.json');

paramobj = BatteryInputParams(jsonstruct);

% Some shortcuts used for the sub-models
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
eac     = 'ElectrodeActiveComponent';
cc      = 'CurrentCollector';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';

%% We setup the battery geometry.
% Here, we use a 1D model and the class BatteryGenerator1D already contains the discretization parameters
gen = BatteryGenerator1D();
gen.fac = 10;
gen = gen.applyResolutionFactors();

% We update pamobj with grid data
paramobj = gen.updateBatteryInputParams(paramobj);

% In this case, we change some of the values of the paramaters that were given in the json file to other values. This is
% done directly on the object paramobj.
paramobj.(ne).(cc).EffectiveElectricalConductivity = 1e5;
paramobj.(pe).(cc).EffectiveElectricalConductivity = 1e5;
paramobj.(thermal).externalTemperature = paramobj.initT;

%%  The Battery model is initialized by sending paramobj to the Battery class constructor 
% see :class:`Battery <Battery.Battery>`

model = Battery(paramobj,'use_thermal',true,'use_solid_diffusion',true);

model.AutoDiffBackend= AutoDiffBackend();

%% We compute the cell capacity and chose a discharge rate
C      = computeCellCapacity(model);
CRate  = 1; 
inputI = (C/hour)*CRate; % current 

%% We setup the schedule 
% We use different time step for the activation phase (small time steps) and the following discharging phase

% We start with rampup time steps to go through the activation phase 
fac=2;
total = 1.4*hour/CRate;
n=10;
dt0=total*1e-6;
times = getTimeSteps(dt0,n, total,fac);
dt= diff(times);
step = struct('val',diff(times),'control',ones(size(dt)));


% We set up a stopping function. Here, the simulation will stop if the output voltage reach a value smaller than 2. This
% stopping function will not be triggered in this case as we switch to voltage control when E=3.6 (see value of inputE
% below).
pe = 'PositiveElectrode';
cc = 'CurrentCollector';
stopFunc = @(model, state, state_prev) (state.(pe).(cc).E < 2.0); 

tup = 0.1; % rampup value for the current function, see rampupSwitchControl
inputE = 3.0; % Value when current control switches to voltage control
srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, inputI, inputE);

% we setup the control by assigning a source and stop function.
control = repmat(struct('src', srcfunc, 'stopFunction', stopFunc), 1, 1); 

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step); 

%%  We setup the initial state
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

%%  We process output and recover the output voltage and current from the output states.
ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
Enew = cellfun(@(x) x.(pe).(cc).E, states); 
Inew = cellfun(@(x) x.(pe).(cc).I, states);
Tmax = cellfun(@(x) max(x.ThermalModel.T), states);
[SOCN,SOCP] =  cellfun(@(x) model.calculateSOC(x), states);
time = cellfun(@(x) x.time, states); 

%% We plot the the output voltage and current

figure
plot((time/hour), Enew, '*-', 'linewidth', 3)
title('Potential (E)')
xlabel('time (hours)')

figure
plot((time/hour), Inew, '*-', 'linewidth', 3)
title('Current (I)')
xlabel('time (hours)')

figure
plot((time/hour), Tmax, '*-', 'linewidth', 3)
title('max(T)')
xlabel('time (hours)')

figure
plot((time/hour), [SOCP,SOCN], '*-', 'linewidth', 3)
title('SOC')
xlabel('time (hours)')
legend('SOC positive','SOC negative')



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
