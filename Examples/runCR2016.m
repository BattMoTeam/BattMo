clear all
close all


% Setup mrst modules

mrstModule add ad-core mrst-gui mpfa agmg



%% Shorthand
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
am      = 'ActiveMaterial';
cc      = 'CurrentCollector';
elyte   = 'Electrolyte';
sep     = 'Separator';
thermal = 'ThermalModel';
ctrl    = 'Control';

%% Setup the geometrical parameters for a CR2016 battery. 
CR2016_diameter = 20*milli*meter;
CR2016_thickness = 1.6*milli*meter;

% Thicknesses of each component
thickness = containers.Map();
% thickness('PositiveCurrentCollector') = 0.1*milli*meter;
% thickness('PositiveActiveMaterial')   = 0.45*milli*meter;
% thickness('ElectrolyteSeparator')     = 0.25*milli*meter;
% thickness('NegativeActiveMaterial')   = 0.7*milli*meter;
% thickness('NegativeCurrentCollector') = thickness('PositiveCurrentCollector');

% Proportions from runBattery3D
zlength = [10; 100; 50; 80; 10];
zz = CR2016_thickness * zlength / sum(zlength);
thickness('PositiveCurrentCollector') = zz(1);
thickness('PositiveActiveMaterial')   = zz(2);
thickness('ElectrolyteSeparator')     = zz(3);
thickness('NegativeActiveMaterial')   = zz(4);
thickness('NegativeCurrentCollector') = thickness('PositiveCurrentCollector');



assert(abs(sum(cell2mat(thickness.values)) - CR2016_thickness) < eps);

% Diameters of each component
diameter = containers.Map();
% diameter('PositiveCurrentCollector') = CR2016_diameter;
% diameter('PositiveActiveMaterial')   = 16*milli*meter;
% diameter('ElectrolyteSeparator')     = 18*milli*meter;
% diameter('NegativeActiveMaterial')   = diameter('PositiveActiveMaterial');
% diameter('NegativeCurrentCollector') = diameter('PositiveCurrentCollector');
diameter('PositiveCurrentCollector') = CR2016_diameter;
diameter('PositiveActiveMaterial')   = diameter('PositiveCurrentCollector');
diameter('ElectrolyteSeparator')     = diameter('PositiveCurrentCollector');
diameter('NegativeActiveMaterial')   = diameter('PositiveCurrentCollector');
diameter('NegativeCurrentCollector') = diameter('PositiveCurrentCollector');

% Angle of the sector
angle = pi / 20;

% Chamfer (not implemented)
chamfer = CR2016_thickness / 6;

% Possible offset (cannot be used for sector grid)
offset = [0, 0]; 

nR = 10;

% Discretization cells in each layer
nLayer = containers.Map();
% nLayer('PositiveCurrentCollector') = 1;
% nLayer('PositiveActiveMaterial')   = 3;
% nLayer('ElectrolyteSeparator')     = 2; 
% nLayer('NegativeActiveMaterial')   = 3;
% nLayer('NegativeCurrentCollector') = nLayer('PositiveCurrentCollector');
hz = min(cell2mat(thickness.values))/2;
for k = keys(thickness)
    key = k{1};
    nLayer(key) = round(thickness(key)/hz);
end

% Discretization cells radially
%offset = [1, 1]*milli*meter;
offset = [];

params = struct('thickness', thickness, ...
                'diameter', diameter, ...
                'nLayer', nLayer, ...
                'offset', offset, ...
                'angle', angle, ...
                'nR', nR);

jsonstruct = parseBattmoJson('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json');
paramobj = BatteryInputParams(jsonstruct); 

gen = CoinCellBatteryGenerator();
%gen = CoinCellSectorBatteryGenerator();
paramobj = gen.updateBatteryInputParams(paramobj, params);

% FIXME
paramobj.(ne).(cc).EffectiveElectricalConductivity = 1e5;
paramobj.(pe).(cc).EffectiveElectricalConductivity = 1e5;

model = Battery(paramobj); 
model.AutoDiffBackend = AutoDiffBackend();




%% Plot the mesh
% The mesh is plotted using the plotGrid() function from MRST. 
colors = crameri('vik', 5);
figure
plotGrid(model.(ne).(cc).G,     'facecolor', colors(1,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
plotGrid(model.(ne).(am).G,     'facecolor', colors(2,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
plotGrid(model.(elyte).(sep).G, 'facecolor', colors(3,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
plotGrid(model.(pe).(am).G,     'facecolor', colors(4,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
plotGrid(model.(pe).(cc).G,     'facecolor', colors(5,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
axis tight;
legend({'negative electrode current collector' , ...
        'negative electrode active material'   , ...
        'separator'                            , ...
        'positive electrode active material'   , ...
        'positive electrode current collector'}, ...
       'location', 'southwest')
view(3)
drawnow


% return


%% C-rate
%CRate = model.(ctrl).CRate;
CRate = 1;

%% Time step schedule
n         = 25; 
dt        = []; 
dt        = [dt; repmat(0.5e-4, n, 1).*1.5.^[1:n]']; 
totalTime = 1.4*hour/CRate;
n         = 40; 
dt        = [dt; repmat(totalTime/n, n, 1)]; 
times     = [0; cumsum(dt)]; 
tt        = times(2 : end); 
step      = struct('val', diff(times), 'control', ones(numel(tt), 1)); 

%% Setup operating limits
tup         = 0.1; % rampup value for the current function, see rampupSwitchControl
srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                            model.Control.Imax, ...
                                            model.Control.lowerCutoffVoltage);
% we setup the control by assigning a source and stop function.
control = struct('src', srcfunc, 'IEswitch', true);

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step); 


%% Initial state
initstate = model.setupInitialState(); 

%% Setup the properties of the nonlinear solver 
nls = NonLinearSolver(); 
% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10; 
nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', {{ctrl,'E'}}, ...
                                                   'targetChangeAbs', 0.03);
% Change default tolerance for nonlinear solver
%model.nonlinearTolerance = 1e-3*model.Control.Imax;
model.nonlinearTolerance = 1e-5; 

nls.errorOnFailure = false; 

% Set verbosity
model.verbose = true;

%% Run the simulation
[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

%%  Process output and recover the output voltage and current from the output states.
ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
Enew = cellfun(@(x) x.(ctrl).E, states); 
Inew = cellfun(@(x) x.(ctrl).I, states);
time = cellfun(@(x) x.time, states); 



%% Plot the mesh
% The mesh is plotted using the plotGrid() function from MRST. 
colors = crameri('vik', 5);
figure
plotGrid(model.(ne).(cc).G,     'facecolor', colors(1,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
plotGrid(model.(ne).(am).G,     'facecolor', colors(2,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
plotGrid(model.(elyte).(sep).G, 'facecolor', colors(3,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
plotGrid(model.(pe).(am).G,     'facecolor', colors(4,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
plotGrid(model.(pe).(cc).G,     'facecolor', colors(5,:), 'edgealpha', 0.5, 'edgecolor', [1, 1, 1]);
axis tight;
legend({'negative electrode current collector' , ...
        'negative electrode active material'   , ...
        'separator'                            , ...
        'positive electrode active material'   , ...
        'positive electrode current collector'}, ...
       'location', 'southwest')
view(3)


figure
plot(time, Inew), title('I')
figure
plot(time, Enew), title('E')

%%
figure
for k = 1:numel(states)
    c{k} = states{k}.PositiveElectrode.ActiveMaterial.c;
end
plotToolbar(model.PositiveElectrode.ActiveMaterial.G, c)
view(3), colorbar
title('PE AM c')

figure
for k = 1:numel(states)
    c{k} = states{k}.PositiveElectrode.ActiveMaterial.phi;
end
plotToolbar(model.PositiveElectrode.ActiveMaterial.G, c)
view(3), colorbar
title('PE AM phi')

figure
for k = 1:numel(states)
    c{k} = states{k}.NegativeElectrode.ActiveMaterial.c;
end
plotToolbar(model.NegativeElectrode.ActiveMaterial.G, c)
view(3), colorbar
title('NE AM c')
figure
for k = 1:numel(states)
    c{k} = states{k}.NegativeElectrode.ActiveMaterial.phi;
end
plotToolbar(model.NegativeElectrode.ActiveMaterial.G, c)
view(3), colorbar
title('NE AM phi')

figure
for k = 1:numel(states)
    c{k} = states{k}.ThermalModel.T;
end
plotToolbar(model.ThermalModel.G, c)
view(3), colorbar
title('T')


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
