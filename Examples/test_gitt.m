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

% Cell capacity
C = computeCellCapacity(model);
% C Rate
CRate = 1;
inputI  = (C/hour)*CRate;
inputE = 4.2;

% First phase : activation phase (with rampup) and charging to inputE
tup = 0.1;
srcfunc_init = @(time, I, E) rampupSwitchControl(time, tup, I, E, inputI, inputE);

% Second phase : Alternate discharging and relaxation 
srcfunc_dis = @(time, I, E) Icontrol(I, inputI);
srcfunc_rel = @(time, I, E) Icontrol(I, 0);

control(1) = struct('src', srcfunc_init, 'stopFunction', stopFunc); 
control(2) = struct('src', srcfunc_dis, 'stopFunction', stopFunc);
control(3) = struct('src', srcfunc_rel, 'stopFunction', stopFunc);

% Setup time steps
n_init = 10;
dt_init = rampupTimesteps(1.5*tup, tup/n_init, 10);
controlInd_init = ones(numel(dt_init), 1);

time_dis = 2*minute;
time_rel = 20*minute;
n_dis = 3;
n_rel = 3;

dt_dis = time_dis/n_dis*ones(n_dis, 1);
dt_rel = time_rel/n_rel;
dt_rel = rampupTimesteps(time_rel, dt_rel, 3);
n_rel = numel(dt_rel);

dt_cycle = [dt_dis; dt_rel];
controlInd_cycle = [2*ones(n_dis, 1); 3*ones(n_rel, 1)];

% N : number of discharging cycles
N = 40;

step.val = [dt_init; ...
            repmat(dt_cycle, N, 1)];
step.control = [controlInd_init; ...
                repmat(controlInd_cycle, N, 1)];

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