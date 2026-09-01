%% Chen model
% Include presentation of the test case (use rst format)

% clear the workspace and close open figures
clear
close all

% load MRST modules
mrstModule add ad-core mrst-gui mpfa

% Useful abbreviations
elyte   = 'Electrolyte';
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
co      = 'Coating';
cc      = 'CurrentCollector';
am      = 'ActiveMaterial';
itf     = 'Interface';
sd      = 'SolidDiffusion';
thermal = 'ThermalModel';

%% model setup
% We load the json input files and compose them

% material properties
jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
jsonstruct_material = removeJsonStructField(jsonstruct_material, {'Control'});

% geometry
jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));

% control
jsonstruct_control = parseBattmoJson('chenexample_control.json');

% We merge the json structures above
jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                               jsonstruct_geometry, ...
                               jsonstruct_control});

% We setup the model but do not run it as we want to change the initial state to compare with PyBamm experiment.
output = runBatteryJson(jsonstruct, 'runSimulation', false);

model    = output.model;
schedule = output.schedule;
nls      = output.nonLinearSolver

%%  We setup the initial state

c_ne = 29.866*mol/litre; % initial concentration at negative electrode
c_pe = 17.038*mol/litre; % initial concentration at positive electrode

% we use a dedicated function. There is no generic function for that
initstate = initStateChen2020(model, c_ne, c_pe);
          
%% We run the simulation

[~, states, ~] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

%% Simplified diffusion model
% We want to consider the simplified diffusion model and therefore setup a model that can run it.

jsonstruct.(ne).(co).(am).diffusionModelType = 'simple';
jsonstruct.(pe).(co).(am).diffusionModelType = 'simple';

output = runBatteryJson(jsonstruct, 'runSimulation', false);

% note that the values for the schedule and non linear solver are unchanged so that we only retrieve the new model for simplified diffusion process.
model2 = output.model;

initstate2 = initStateChen2020(model2, c_ne, c_pe);

% Run simulation
[~, states2, ~] = simulateScheduleAD(initstate2, model2, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

%%  We process output and recover the output voltage and current from the output states.

ind     = cellfun(@(x) not(isempty(x)), states);
states  = states(ind);
E1      = cellfun(@(state) state.Control.E, states);
I1      = cellfun(@(state) state.Control.I, states);
time1   = cellfun(@(state) state.time, states);
time1   = time1/hour;

ind     = cellfun(@(x) not(isempty(x)), states2);
states2 = states2(ind);
E2      = cellfun(@(state) state.Control.E, states2);
I2      = cellfun(@(state) state.Control.I, states2);
time2   = cellfun(@(state) state.time, states2);
time2   = time2/hour;

%% We plot the the output voltage and current and compare with PyBaMM results

% We load the PyBaMM results
loadChenPybammSolution
[t1, u1] = deal(t, u);
[t2, u2] = deal(t_infdiff, u_infdiff);

doplot = true;

if doplot

    % We plot the solutions
    l = lines(4);
    figure
    hold on
    plot(t1, u1, 'linewidth', 3, 'color', l(1, :), 'linestyle', '--', 'displayname', 'pybamm - solid diffusion')
    plot(t2, u2, 'linewidth', 3, 'color', l(2, :), 'linestyle', '--', 'displayname', 'pybamm - instantaneous solid diffusion')
    plot(time1, E1,'-', 'linewidth', 3, 'color', l(3, :), 'displayname', 'battmo')
    plot(time2, E2,'-', 'linewidth', 3, 'color', l(4, :), 'displayname', 'battmo - simplified diffusion model')

    set(gca, 'fontsize', 18);
    title('Cell Voltage / V')
    xlabel('time (hours)')
    legend('fontsize', 18, 'location', 'southwest')

end

% Make sure BattMo and PyBaMM match
N = 1000;
x = linspace(time1(1), t1(end), N);
pb = interp1(t1, u1, x);
battmo = interp1(time1, E1, x);
assert(norm(pb - battmo) / norm(pb) < 1e-3);



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
