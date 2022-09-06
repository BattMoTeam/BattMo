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
jsonstruct = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));

paramobj = BatteryInputParams(jsonstruct);


%% We setup the battery geometry ("bare" battery with no current collector).

gen = BareBatteryGenerator3D();

% We update paramobj with grid data
paramobj = gen.updateBatteryInputParams(paramobj);


%%  The Battery model is initialized by sending paramobj to the Battery class constructor 

model = Battery(paramobj);

%% We fix the input current to 5A

model.Control.Imax = 5*ampere;

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

tup = 0.1; % rampup value for the current function, see rampupSwitchControl
srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                            model.Control.Imax, ...
                                            model.Control.lowerCutoffVoltage);
% we setup the control by assigning a source and stop function.
control = struct('src', srcfunc, 'IEswitch', true);

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step); 

%%  We setup the initial state

c_ne = 29.866*mol/litre; % initial concentration at negative electrode
c_pe = 17.038*mol/litre; % initial concentration at positive electrode
initstate = initStateChen2020(model, c_ne, c_pe);

%% We setup the solver
% Default values would typically work.

nls = NonLinearSolver(); 
% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10; 
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false;
linearsolver = 'direct';
switch linearsolver
  case 'agmg'
    mrstModule add agmg
    nls.LinearSolver = AGMGSolverAD('verbose', true, 'reduceToCell', true); 
    nls.LinearSolver.tolerance = 1e-3; 
    nls.LinearSolver.maxIterations = 30; 
    nls.maxIterations = 10; 
    nls.verbose = 10;
  case 'battery'
    nls.LinearSolver = LinearSolverBatteryExtra('verbose', false, 'reduceToCell', false, 'verbosity', 3, 'reuse_setup', false, 'method', 'direct');
    nls.LinearSolver.tolerance=0.5e-4*2;          
  case 'direct'
    % nothing to setup in this case
  otherwise
    error('linear solver not recognized.');
end

% Change default tolerance for nonlinear solver
model.nonlinearTolerance = 1e-5; 
% Set verbosity
model.verbose = false;

model.AutoDiffBackend= AutoDiffBackend();

%% We run simulation

[~, states, ~] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

% We want to consider the simplified diffusion model and therefore setup a model that can run it.

jsonstruct.NegativeElectrode.ActiveMaterial.diffusionModelType = 'simple';
jsonstruct.PositiveElectrode.ActiveMaterial.diffusionModelType = 'simple';

paramobj = BatteryInputParams(jsonstruct);
paramobj = gen.updateBatteryInputParams(paramobj);

model2 = Battery(paramobj);

% setup initial state (the variables are different for the simplified model and therfore the initialization is not the same)

initstate2 = initStateChen2020(model2, c_ne, c_pe);

% Run simulation
[~, states2, ~] = simulateScheduleAD(initstate2, model2, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);


%%  We process output and recover the output voltage and current from the output states.

E1    = cellfun(@(state) state.Control.E, states); 
I1    = cellfun(@(state) state.Control.I, states);
time1 = cellfun(@(state) state.time, states);

E2    = cellfun(@(state) state.Control.E, states2); 
I2    = cellfun(@(state) state.Control.I, states2);
time2 = cellfun(@(state) state.time, states2);


%% We plot the the output voltage and current and compare with pybamm results

% We load the pybamm results
loadChenPybammSolution
[t1, u1] = deal(t, u);
[t2, u2] = deal(t_infdiff, u_infdiff);

% We plot the solutions

set(0, 'defaultFigurePosition', [671 510 900 600]);

l = lines(4);
figure
hold on
plot(t1, u1, 'linewidth', 3, 'color', l(1, :), 'displayname', 'pybamm - solid diffusion')
plot(t2, u2, 'linewidth', 3, 'color', l(2, :), 'displayname', 'pybamm - instantaneous solid diffusion')
plot((time1/hour), E1,'-', 'linewidth', 3, 'color', l(3, :), 'displayname', 'battmo')
plot((time2/hour), E2,'-', 'linewidth', 3, 'color', l(4, :), 'displayname', 'battmo - simplified diffusion model')


set(gca, 'fontsize', 18);
title('Cell Voltage / V')
xlabel('time (hours)')
legend('fontsize', 18, 'location', 'south west')


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
