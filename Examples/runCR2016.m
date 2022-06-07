clear all
close all


% Setup mrst modules

mrstModule add ad-core mrst-gui mpfa agmg

%% We setup the geometrical parameters for a CR2016 battery. 

% Thickness of each component ordered as
% - positive current collector
% - positive electrode
% - electrolyte separator 
% - negative electrode
% - negative current collector

CR2016_diameter = 20*milli*meter;
CR2016_thickness = 1.6*milli*meter;

% Thicknesses of each component
thickness = containers.Map();
thickness('PositiveCurrentCollector') = 0.1*milli*meter;
thickness('PositiveActiveMaterial')   = 0.7*milli*meter;
thickness('ElectrolyteSeparator')     = 0.25*milli*meter;
thickness('NegativeActiveMaterial')   = 0.45*milli*meter;
thickness('NegativeCurrentCollector') = thickness('PositiveCurrentCollector');
assert(abs(sum(cell2mat(thickness.values)) - CR2016_thickness) < eps);

% Diameters of each component
diameter = containers.Map();
diameter('PositiveCurrentCollector') = CR2016_diameter;
diameter('PositiveActiveMaterial')   = 16*milli*meter;
diameter('ElectrolyteSeparator')     = 18*milli*meter;
diameter('NegativeActiveMaterial')   = diameter('PositiveActiveMaterial');
diameter('NegativeCurrentCollector') = diameter('PositiveCurrentCollector');

% Angle of the sector
angle = pi / 20;

% % Chamfer (not implemented)
% chamfer = CR2016_thickness / 6;

% % Possible offset (cannot be used for sector grid)
% offset = [0, 0]; 

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
nR = 10;

params = struct('thickness', thickness, ...
                'diameter', diameter, ...
                'angle', angle, ...
                'nLayer', nLayer, ...
                'nR', nR);

jsonstruct = parseBattmoJson('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json');
paramobj = BatteryInputParams(jsonstruct); 

gen = CoinCellSectorBatteryGenerator();
paramobj = gen.updateBatteryInputParams(paramobj, params);

model = Battery(paramobj); 
model.AutoDiffBackend = AutoDiffBackend();


%% Shorthand
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';

%% C-rate
CRate = model.(ctrl).CRate;

%% Time step schedule
switch model.(ctrl).controlPolicy
  case 'CCCV'
    total = 3.5*hour/CRate;
  case 'IEswitch'
    total = 1.4*hour/CRate;
  otherwise
    error('control policy not recognized');
end

n     = 100;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

switch model.Control.controlPolicy
  case 'IEswitch'
    tup = 0.1; % rampup value for the current function, see rampupSwitchControl
    srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                model.Control.Imax, ...
                                                model.Control.lowerCutoffVoltage);
    % we setup the control by assigning a source and stop function.
    control = struct('src', srcfunc, 'IEswitch', true);
  case 'CCCV'
    control = struct('CCCV', true);
  otherwise
    error('control policy not recognized');
end

schedule = struct('control', control, 'step', step); 

%% Initial state
initstate = model.setupInitialState(); 

%% Setup the properties of the nonlinear solver 
nls = NonLinearSolver(); 
% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10; 
nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);
% Change default tolerance for nonlinear solver
model.nonlinearTolerance = 1e-3*model.Control.Imax;
% Set verbosity
model.verbose = true;

%% Run the simulation
[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 




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
