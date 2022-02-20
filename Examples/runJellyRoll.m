% Setup mrst modules

mrstModule add ad-core multimodel mrst-gui mpfa agmg

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
nas = 10; 

% Number of discretization cells in the longitudonal
nL = 3;

% structure that describes the tab setups (see SpiralBatteryGenerator)
tabparams.tabcase   = 'aligned tabs';
tabparams.width     = 3*milli*meter;
tabparams.fractions = linspace(0.01, 0.9, 6);

spiralparams = struct('nwindings'   , nwindings, ...
                      'rInner'      , rInner   , ...
                      'widthDict'   , widthDict, ...
                      'nrDict'      , nrDict   , ...
                      'nas'         , nas      , ...
                      'L'           , L        , ...
                      'nL'          , nL       , ...
                      'tabparams'   , tabparams, ...
                      'angleuniform', false); 

% The input material parameters given in json format are used to populate the paramobj object.
jsonfilename ='JsonDatas/lithiumbattery.json';
jsonstruct = parseBatmoJson(jsonfilename); 
paramobj = BatteryInputParams(jsonstruct); 

th = 'ThermalModel';
paramobj.(th).externalHeatTransferCoefficientSideFaces = 100*watt/meter^2;
paramobj.(th).externalHeatTransferCoefficientTopFaces = 10*watt/meter^2;

gen = SpiralBatteryGenerator(); 

paramobj = gen.updateBatteryInputParams(paramobj, spiralparams);

model = Battery(paramobj); 

[cap, cap_neg, cap_pos, specificEnergy] = computeCellCapacity(model);
fprintf('ratio : %g, energy : %g\n', cap_neg/cap_pos, specificEnergy/hour);

%% Setup schedule
CRate = 1;

fac   = 2; 
total = 1.4*hour/CRate; 
n     = 10; 
dt0   = total*1e-6; 
times = getTimeSteps(dt0, n, total, fac); 

%% We compute the cell capacity, which used to compute schedule from CRate
C = computeCellCapacity(model); 
inputI = (C/hour)*CRate; 
inputE = 3; 

tt = times(2 : end); 

step = struct('val', diff(times), 'control', ones(numel(tt), 1)); 

pe = 'PositiveElectrode'; 
cc = 'CurrentCollector'; 
stopFunc = @(model, state, state_prev) (state.(pe).(cc).I < 1e-3*inputI && state.time> hour/(2*CRate)); 
tup = 0.1/CRate; 

simcase = 'discharge';

switch simcase
    
  case 'discharge'
    stopFunc = @(model, state, state_prev) (state.(pe).(cc).E < inputE+1e-4); 
    srcfunc   = @(time, I, E) rampupSwitchControl(time, tup, I, E, inputI, inputE); 
    control   = repmat(struct('src', srcfunc, 'stopFunction', stopFunc), 1, 1); 
    schedule  = struct('control', control, 'step', step); 

    %% We setup the initial state
    initstate = model.setupInitialState(); 
    initElytec = 1*mol/litre;
    initstate.Electrolyte.c = initElytec*ones(model.Electrolyte.G.cells.num, 1);
    
  case 'charge'
    
    model.SOC = 0.01;
    initstate = model.setupInitialState();
    stopFunc = @(model, state, state_prev) (state.(pe).(cc).I > - 1e-3*inputI  && state.time> hour/(2*CRate)); 
    srcfunc  = @(time, I, E) rampupSwitchControl(time, tup, I, E, -inputI, 4.2); 
    control  = repmat(struct('src', srcfunc, 'stopFunction', stopFunc), 1, 1); 
    schedule = struct('control', control, 'step', step); 
    
  otherwise
    error('simcase not recognized')

end

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

nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', {{'PositiveElectrode', 'CurrentCollector', 'E'}}, 'targetChangeAbs', 0.03);
linearsolver = 'direct';
switch linearsolver
  case 'agmg'
    mrstModule add agmg
    nls.LinearSolver = AGMGSolverAD('verbose', true, 'reduceToCell', true); 
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
problem = packSimulationProblem(initstate, model, schedule, dataFolder, 'Name', 'jellyroll', 'NonLinearSolver', nls);
problem.SimulatorSetup.OutputMinisteps = true; 

resetSimulation = false;
if resetSimulation
    %% clear previously computed simulation
    clearPackedSimulatorOutput(problem);
end
simulatePackedProblem(problem);
[globvars, states, report] = getPackedSimulatorOutput(problem);






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