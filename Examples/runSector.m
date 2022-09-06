% Setup mrst modules

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


% Number of discretization cells in the longitudonal
nL  = 5;
nas = 3;
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
                      'tabparams'   , [], ...
                      'angleuniform', true); 

% The input material parameters given in json format are used to populate the paramobj object.
jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
jsonstruct.include_current_collectors = true;

paramobj = BatteryInputParams(jsonstruct); 

th = 'ThermalModel';
%paramobj.(th).externalHeatTransferCoefficientSideFaces = 100*watt/meter^2;
%paramobj.(th).externalHeatTransferCoefficientTopFaces = 10*watt/meter^2;
paramobj.(th).externalHeatTransferCoefficient=10;
%gen = SpiralBatteryGenerator(); 
gen = SectorBatteryGenerator(); 

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

tup = 0.1/CRate; 

simcase = 'discharge';

switch simcase
    
  case 'discharge'
    srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                model.Control.Imax, ...
                                                model.Control.lowerCutoffVoltage);
    % we setup the control by assigning a source and stop function.
    control = struct('src', srcfunc, 'IEswitch', true);
    schedule  = struct('control', control, 'step', step); 

    %% We setup the initial state
    initstate = model.setupInitialState(); 
    initElytec = 1*mol/litre;
    initstate.Electrolyte.c = initElytec*ones(model.Electrolyte.G.cells.num, 1);
    
  case 'charge'
    
    model.SOC = 0.01;
    initstate = model.setupInitialState();
    srcfunc  = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                 - model.Control.Imax, ...
                                                 model.Control.lowerCutoffVoltage); 
    control = struct('src', srcfunc, 'IEswitch', true);
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

nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', {{'Control', 'E'}}, 'targetChangeAbs', 0.03);
linearsolver = 'direct';
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
dataFolder = 'BattMo';
problem = packSimulationProblem(initstate, model, schedule, dataFolder, 'Name', 'jellyroll', 'NonLinearSolver', nls);
problem.SimulatorSetup.OutputMinisteps = true; 

resetSimulation = true;
if resetSimulation
    %% clear previously computed simulation
    clearPackedSimulatorOutput(problem, 'prompt', false);
end
simulatePackedProblem(problem);
[globvars, states, report] = getPackedSimulatorOutput(problem);


%% plot thermal

figure
plotToolbar(model.G, states);
view([0,-1,0]);

%%  Process output and recover the output voltage and current from the output states.

ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
ctrl = 'Control';
E = cellfun(@(x) x.(ctrl).E, states); 
I = cellfun(@(x) x.(ctrl).I, states);
time = cellfun(@(x) x.time, states); 

figure
plot(time, E);

%% Plot an animated summary of the results
%plotDashboard(model, states, 'step', 0);



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
