%% Battery 1D model
% Include presentation of the test case (use rst format)
clear all

% load MRST modules
mrstModule add ad-core multimodel mrst-gui mpfa

% We create an instance of BatteryInputParams. This class is used to initiate the battery simulator and it propagates
% all the parameters through out the submodels.

% The input parameters can be given in json format. The json file is read and used to populate the paramobj object.
jsonstruct = parseBatmoJson('JsonDatas/Chen2020/chenBattery.json');

paramobj = BareBatteryInputParams(jsonstruct);

amn = paramobj.NegativeElectrode.ActiveMaterial;
amp = paramobj.PositiveElectrode.ActiveMaterial;

amn.volumetricSurfaceArea = 1;
amp.volumetricSurfaceArea = 1;

amn.rp = (amn.rp^2)/(3*amn.volumeFraction);
amp.rp = (amp.rp^2)/(3*amp.volumeFraction);

paramobj.NegativeElectrode.ActiveMaterial = amn;
paramobj.PositiveElectrode.ActiveMaterial = amp;

% Some shortcuts used for the sub-models
ne    = 'NegativeElectrode';
pe    = 'PositiveElectrode';
elyte = 'Electrolyte';

%% We setup the battery geometry.
% Here, we use a 1D model and the class BatteryGenerator1D already contains the discretization parameters
gen = BareBatteryGenerator3D();
% We update pamobj with grid data
paramobj = gen.updateBatteryInputParams(paramobj);

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
dt1   = rampupTimesteps(0.1, 0.1, 20);
% We choose time steps for the rest of the simulation (discharge phase)
dt2   = 2e-2*hour*ones(200, 1);
% We concatenate the time steps
dt    = [dt1; dt2];
times = [0; cumsum(dt)]; 
tt    = times(2 : end); 
step  = struct('val', diff(times), 'control', ones(numel(tt), 1)); 

stopFunc = @(model, state, state_prev) (state.(pe).E < 2.0); 

tup = 0.1; % rampup value for the current function, see rampupSwitchControl
inputE = 3.6; % Value when current control switches to voltage control
srcfunc = @(time) rampupControl(time, tup, inputI);

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
cs = cell(2,1);
initstate.(elyte).cs = cs;
initstate.(elyte).cs{1} = 1000*ones(bat.(elyte).G.cells.num, 1);

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
model.verbose = false;

% Run simulation

[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 


%%  We process output and recover the output voltage and current from the output states.
ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
Enew = cellfun(@(x) x.(pe).E, states); 
Inew = cellfun(@(x) x.(pe).I, states);
time = cellfun(@(x) x.time, states); 

%% We plot the the output voltage and current

figure
plot((time/hour), Enew, '*-', 'linewidth', 3)
title('Potential (E)')
xlabel('time (hours)')
hold on
loadChenPybammSolution
plot(t, u, 'linewidth', 3, 'color', 'green')
legend({'batmo', 'pybamm'})
figure
plot((time/hour), Inew, '*-', 'linewidth', 3)
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
