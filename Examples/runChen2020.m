%% Chen model
% Include presentation of the test case (use rst format)

% clear the workspace and close open figures
clear
close all
clc

% load MRST modules
mrstModule add ad-core mrst-gui mpfa

% We create an instance of BatteryInputParams. This class is used to initiate the battery simulator and it propagates
% all the parameters through out the submodels.

% The input parameters can be given in json format. The json file is read and used to populate the paramobj object.
jsonstruct = parseBatmoJson('JsonDatas/Chen2020/chenBattery.json');

paramobj = BareBatteryInputParams(jsonstruct);

% Some shortcuts used for the sub-models
ne    = 'NegativeElectrode';
pe    = 'PositiveElectrode';
elyte = 'Electrolyte';

%% We setup the battery geometry.
gen = BareBatteryGenerator3D();
% We update pamobj with grid data
paramobj = gen.updateBatteryInputParams(paramobj);

paramobj.NegativeElectrode.InterDiffusionCoefficient = 0;
paramobj.PositiveElectrode.InterDiffusionCoefficient = 0;

%%  The Battery model is initialized by sending paramobj to the Battery class constructor 

model = BareBattery(paramobj);

%% We compute the cell capacity and chose a discharge rate
% C      = computeCellCapacity(model);
% CRate  = 1/5; 
% inputI = (C/hour)*CRate; % current 
inputI = 5;

%% We setup the schedule 
% We use different time step for the activation phase (small time steps) and the following discharging phase
% We start with rampup time steps to go through the activation phase 

fac   = 2; 
total = 1.4*hour; 
n     = 100; 
dt0   = total*1e-6; 
times = getTimeSteps(dt0, n, total, fac); 
dt    = diff(times);
dt    = dt(1 : end);
step  = struct('val', dt, 'control', ones(size(dt)));

% We set up a stopping function. Here, the simulation will stop if the output voltage reach a value smaller than 2. This
% stopping function will not be triggered in this case as we switch to voltage control when E=3.6 (see value of inputE
% below).
pe = 'PositiveElectrode';
cc = 'CurrentCollector';
stopFunc = @(model, state, state_prev) (state.(pe).E < 2.6); 

tup     = 0.1; % rampup value for the current function, see rampupSwitchControl
inputE  = 2;   % Value when current control switches to voltage control
srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, inputI, inputE);

% we setup the control by assigning a source and stop function.
control = repmat(struct('src', srcfunc, 'stopFunction', stopFunc), 1, 1); 

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step); 


%%  We setup the initial state

nc = model.G.cells.num;
T = model.initT;
initstate.T = T*ones(nc, 1);

bat = model;
elyte = 'Electrolyte';
ne    = 'NegativeElectrode';
pe    = 'PositiveElectrode';
am    = 'ActiveMaterial';

initstate = model.updateTemperature(initstate);

% we setup negative electrode initial state
nam = bat.(ne).(am); 

c = 29866.0;
c = c*ones(nam.G.cells.num, 1);
initstate.(ne).c = c;
% We bypass the solid diffusion equation to set directly the particle surface concentration (this is a bit hacky)
initstate.(ne).(am).cElectrode = c;
initstate.(ne).(am) = nam.updateOCP(initstate.(ne).(am));

OCP = initstate.(ne).(am).OCP;
ref = OCP(1);

initstate.(ne).phi = OCP - ref;

% we setup positive electrode initial state
pam = bat.(pe).(am); 

c = 17038.0;
c = c*ones(pam.G.cells.num, 1);
initstate.(pe).c = c;
% We bypass the solid diffusion equation to set directly the particle surface concentration (this is a bit hacky)
initstate.(pe).(am).cElectrode = c;
initstate.(pe).(am) = pam.updateOCP(initstate.(pe).(am));

OCP = initstate.(pe).(am).OCP;

initstate.(pe).phi = OCP - ref;

initstate.(elyte).phi = zeros(bat.(elyte).G.cells.num, 1) - ref;
initstate.(elyte).c = 1000*ones(bat.(elyte).G.cells.num, 1);

% setup initial positive electrode external coupling values

initstate.(pe).E = OCP(1) - ref;
initstate.(pe).I = 0;
            
% Setup nonlinear solver 
nls = NonLinearSolver(); 
% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10; 
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false; 
% Change default tolerance for nonlinear solver
model.nonlinearTolerance = 1e-5; 
% Set verbosity
model.verbose = true;

% Run simulation
[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

%%  We process output and recover the output voltage and current from the output states.
ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
Enew = cellfun(@(x) x.(pe).E, states); 
Inew = cellfun(@(x) x.(pe).I, states);
time = cellfun(@(x) x.time, states); 

%% We plot the the output voltage and current

loadChenPybammSolution
[t1, u1] = deal(t, u);
loadChenPybammSolution2
[t2, u2] = deal(t, u);

run('/home/xavier/Python/PyBaMM/readme.m')

set(0, 'defaultFigurePosition', [671 510 900 600]);

l = lines(3);
figure
plot((time/hour), Enew,'-', 'linewidth', 3, 'color', l(1, :), 'displayname', 'batmo - instantaneous diffusion')
hold on
plot(t1, u1, 'linewidth', 3, 'color', l(2, :), 'displayname', 'pybamm')
plot(t2, u2, 'linewidth', 3, 'color', l(3, :), 'displayname', 'pybamm - instantaneous diffusion')
% plot(t_pybamm/hour, u_pybamm, '-', 'linewidth', 3, 'color', l(3, :), 'displayname', 'pybamm hello')
set(gca, 'fontsize', 18);
title('Potential (E)')
xlabel('time (hours)')
legend('fontsize', 18, 'location', 'south west')

figure
plot((time/hour), Inew, 'linewidth', 3)
title('Current (I)')
xlabel('time (hours)')


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
