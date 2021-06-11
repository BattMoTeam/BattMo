clear all
close all

% setup mrst modules
mrstModule add ad-core multimodel mrst-gui battery mpfa

mrstVerbose off

paramobj = LithiumBatteryInputParams();

% Setup battery
modelcase = '2D';

switch modelcase

  case '1D'

    gen = BatteryGenerator1D();
    paramobj = gen.updateBatteryInputParams(paramobj);
    paramobj.ne.cc.EffectiveElectricalConductivity = 100;
    paramobj.pe.cc.EffectiveElectricalConductivity = 100;
    schedulecase = 3;
    
    paramobj.thermal.externalHeatTransferCoefficient = 1000;
    paramobj.thermal.externalTemperature = paramobj.initT;

  case '2D'

    gen = BatteryGenerator2D();
    paramobj = gen.updateBatteryInputParams(paramobj);
    schedulecase = 1;

    paramobj.ne.cc.EffectiveElectricalConductivity = 1e5;
    paramobj.pe.cc.EffectiveElectricalConductivity = 1e5;
    
    paramobj.thermal.externalTemperature = paramobj.initT;
    paramobj.SOC = 0.99;
    
    tfac = 1; % used in schedule setup
  
  case '3D'

    gen = BatteryGenerator3D();
    
    fac = 1; 
    gen.facx = fac; 
    gen.facy = fac; 
    gen.facz = fac;
    gen = gen.applyResolutionFactors();
    paramobj = gen.updateBatteryInputParams(paramobj);
    
    schedulecase = 5;
    
    paramobj.thermal.externalTemperature = paramobj.initT;
    
end

model = Battery(paramobj);

initstate = model.setupInitialState(); 

pe = 'PositiveElectrode';
cc = 'CurrentCollector';
stopFunc = @(model, state, state_prev) (state.(pe).(cc).E < 2.0); 


%% Gitt inputs

% C Rate
CRate = 2;

% pulseFraction
pulsefraction = 0.01;

% time relaxation
time_rel = 4*hour;

% Rampup time (to cope with the abrupt change on control)
time_rampup = 1*second; % rampup time (linear interpolation between the two states)

% number of cycle (default could be 1/pulfraction)
N = 2;

% Discretization parameters
N_per_dis    = 5; % Number of time step in discharge phase
N_per_rel    = 5; % Number of time step in relaxation phase
N_per_rampup = 3; % Number of time step in rampup phase

%% setup 

% discharge time
time_dis    = pulsefraction*CRate*hour; % time of discharging

% Cell capacity
C = computeCellCapacity(model);
inputI = (C/hour)*CRate;
inputE = 4.2;

% Activtation phase : activation phase (with rampup) and charging to inputE

tup = 0.1; % rampup time
time_init = 5*minute;
srcfunc_init = @(time, I, E) rampupSwitchControl(time, tup, I, E, inputI, inputE);

% Gitt phase : Alternate discharging and relaxation 

% we setup tabulated input for the gitt part

dtpoints = [time_rampup; ...
            time_dis - time_rampup;
            time_rampup;
            time_rel - time_rampup];
dIpoints = [1; 0; -1; 0];

tpoints = [0; cumsum(repmat(dtpoints, N, 1))];
Ipoints = [0; cumsum(repmat(dIpoints, N, 1))];

tpoints = time_init + tpoints; % starts at end of activation phase
Ipoints = inputI*Ipoints; % scales with Iinput

srcfunc_gitt = @(time, I, E) tabulatedIControl(time, tpoints, Ipoints);

control(1) = struct('src', srcfunc_init, 'stopFunction', stopFunc); 
control(2) = struct('src', srcfunc_gitt, 'stopFunction', stopFunc);

n_init = 5;
dt_init = rampupTimesteps(time_init, time_init/n_init, 3);


dt_cycle = [time_rampup/N_per_rampup*ones(N_per_rampup, 1); ...
            (time_dis - time_rampup)/N_per_dis*ones(N_per_dis, 1); ...
            time_rampup/N_per_rampup*ones(N_per_rampup, 1); ...
            (time_rel - time_rampup)/N_per_rel*ones(N_per_rel, 1)];

dt_cycle = repmat(dt_cycle, N, 1);

step.val = [dt_init; dt_cycle];
step.control = [ones(numel(dt_init), 1); ...
                2*ones(numel(dt_cycle), 1)];

schedule = struct('control', control, 'step', step); 

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
end

use_iterative = false; 
if(use_iterative)
    % nls.LinearSolver = LinearSolverBattery('method', 'iterative'); 
    % nls.LinearSolver = LinearSolverBattery('method', 'direct'); 
    mrstModule add agmg
    nls.LinearSolver = LinearSolverBattery('method', 'agmg', 'verbosity', 1);
    nls.LinearSolver.tol = 1e-3;
    nls.verbose = 10
end
model.nonlinearTolerance = 1e-5; 
model.verbose = false;

% Run simulation

doprofiling = false;
if doprofiling
    profile off
    profile on
end

[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule,...
                                                'OutputMinisteps', true,...
                                                'NonLinearSolver', nls); 
if doprofiling
    profile off
    profile report
end 


%%  Process output

ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
Enew = cellfun(@(x) x.(pe).(cc).E, states); 
Inew = cellfun(@(x) x.(pe).(cc).I, states);
time = cellfun(@(x) x.time, states); 

%%

figure
plot((time/hour), Enew, '*-', 'linewidth', 1)
title('Potential (E)')
xlabel('time (hours)')

figure
plot((time/hour), Inew, '*-', 'linewidth', 1)
title('Current (I)')
xlabel('time (hours)')

return

%% more plotting (case dependent)

switch modelcase
  case '1D'
    % plot1D;
  case '2D'
    plotThermal(model, states);
    plot2Dconc;
  case '3D'
    plot3D;
    plot3Dconc;
    plot3Dphi;
end