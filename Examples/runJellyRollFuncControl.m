function problem = runJellyRollFuncNew(options)
% Setup mrst modules
CRate = options.CRate;
simcase = options.simcase;
nas = options.nas;
nL = options.nL;
tabcase = options.tabcase;
mrstModule add ad-core mrst-gui mpfa agmg

%% We setup the geometrical parameters for a 4680 battery. 
%% Those will be gathered in structure spiralparams (see below) and used by SpiralBatteryGenerator to generate the spiral layered geometry of the jelly roll

% Inner radius of the jelly roll
rInner = 2*milli*meter; 

% widths of each component ordered as
% - positive current collector
% - positive electrode
% - electrolyte separator 
% - negative electrode
% - negative current collector

widths = [25, 64, 15, 57, 15]*micro*meter; 

widthDict = containers.Map( ...
    {'ElectrolyteSeparator',... 
     'NegativeActiveMaterial',...
     'NegativeCurrentCollector',...
     'PositiveActiveMaterial',...
     'PositiveCurrentCollector'},...
    widths); 

nwidths = [widthDict('PositiveActiveMaterial');...
           widthDict('PositiveCurrentCollector');...
           widthDict('PositiveActiveMaterial');...
           widthDict('ElectrolyteSeparator');...
           widthDict('NegativeActiveMaterial');...
           widthDict('NegativeCurrentCollector');...
           widthDict('NegativeActiveMaterial');...
           widthDict('ElectrolyteSeparator')]; 

dr = sum(nwidths);

% Outer radius of the jelly roll
rOuter = 46*milli*meter/2;
% Height of the jelly roll
L = 80*milli*meter; 

dR = rOuter - rInner; 
% Computed number of windings
nwindings = ceil(dR/dr);

% number of discretization cells in radial direction for each component.
nrDict = containers.Map( ...
    {'ElectrolyteSeparator',... 
     'NegativeActiveMaterial',...
     'NegativeCurrentCollector',...
     'PositiveActiveMaterial',...
     'PositiveCurrentCollector'},...
    [3, 3, 3, 3, 3]); 

% Number of discretization cells in the angular direction
%nas = 10; 

% Number of discretization cells in the longitudonal
%nL = 3;

% structure that describes the tab setups (see SpiralBatteryGenerator)
tabparams.tabcase   = tabcase;
tabparams.width     = 3*milli*meter;
tabparams.fractions = linspace(0.01, 0.9, 4);

spiralparams = struct('nwindings'   , nwindings, ...
                      'rInner'      , rInner   , ...
                      'widthDict'   , widthDict, ...
                      'nrDict'      , nrDict   , ...
                      'nas'         , nas      , ...
                      'L'           , L        , ...
                      'nL'          , nL       , ...
                      'tabparams'   , tabparams, ...
                      'angleuniform', true); 

% The input material parameters given in json format are used to populate the paramobj object.
jsonstruct = parseBattmoJson('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json');
jsonstruct.Control.controlPolicy = simcase;
jsonstruct.Control.lowerCutoffVoltage=3.0
paramobj = BatteryInputParams(jsonstruct);
paramobj.SOC = 0.02;
 

th = 'ThermalModel';
paramobj.(th).externalHeatTransferCoefficientSideFaces = 100*watt/meter^2;
paramobj.(th).externalHeatTransferCoefficientTopFaces = 10*watt/meter^2;

gen = SpiralBatteryGenerator(); 

paramobj = gen.updateBatteryInputParams(paramobj, spiralparams);

model = Battery(paramobj); 

[cap, cap_neg, cap_pos, specificEnergy] = computeCellCapacity(model);
fprintf('ratio : %g, energy : %g\n', cap_neg/cap_pos, specificEnergy/hour);

%% Setup schedule
%CRate = 1;
CRate = model.Control.CRate;

%% Setup the time step schedule 
% Smaller time steps are used to ramp up the current from zero to its
% operational value. Larger time steps are then used for the normal
% operation.
switch model.Control.controlPolicy
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

% we setup the control by assigning a source and stop function.
% control = struct('CCCV', true); 
%  !!! Change this to an entry in the JSON with better variable names !!!

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

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step); 
initstate = model.setupInitialState();
initstate.Control.I = model.Control.Imax;
% Setup nonlinear solver 
nls = NonLinearSolver(); 

% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10; 
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false; 
% Change default tolerance for nonlinear solver
model.nonlinearTolerance = 1e-4; 

use_diagonal_ad = false;
if(use_diagonal_ad)
    model.AutoDiffBackend = DiagonalAutoDiffBackend(); 
    model.AutoDiffBackend.useMex = true; 
    model.AutoDiffBackend.modifyOperators = true; 
    model.AutoDiffBackend.rowMajor = true; 
    model.AutoDiffBackend.deferredAssembly = false; % error with true for now
else
    model.AutoDiffBackend = AutoDiffBackend(); 
end
nls.timeStepSelector=StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);
linearsolver = 'agmg';
switch linearsolver
  case 'agmg'
    mrstModule add agmg
    nls.LinearSolver = AGMGSolverAD('verbose', false, 'reduceToCell', true); 
    nls.LinearSolver.tolerance = 1e-3; 
    nls.LinearSolver.maxIterations = 30; 
    nls.maxIterations = 10; 
    nls.verbose = 10;
  case 'direct'
    disp('standard direct solver')
  otherwise
    error()
end

model.nonlinearTolerance = 1e-4; 
model.verbose = true; 

% Run simulation
dataFolder = 'BatMo';
filename = ['jellroll_C_',num2str(CRate),'_',simcase,'_nas_',num2str(nas),'_nL_',num2str(nL),'_',tabcase];
problem = packSimulationProblem(initstate, model, schedule, dataFolder, 'Name', filename, 'NonLinearSolver', nls);
problem.SimulatorSetup.OutputMinisteps = true; 
end
%resetSimulation = true;
%if resetSimulation
    %% clear previously computed simulation
%    clearPackedSimulatorOutput(problem);
%end

%simulatePackedProblem(problem);
%[globvars, states, report] = getPackedSimulatorOutput(problem);






%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
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
