clear all
close all

% setup mrst modules
mrstModule add ad-core multimodel mrst-gui battery mpfa

mrstVerbose off

% Value used in rampup function, see currentSource.
tup = 0.1;

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
    
    paramobj.J = 1e1;
    paramobj.thermal.externalHeatTransferCoefficient = 1000;
    paramobj.thermal.externalTemperature = paramobj.initT;

  case '2D'

    gen = BatteryGenerator2D();
    paramobj = gen.updateBatteryInputParams(paramobj);
    schedulecase = 1;

    paramobj.ne.cc.EffectiveElectricalConductivity = 1e6;
    paramobj.pe.cc.EffectiveElectricalConductivity = 1e6;
    
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
    paramobj.J = 5e-4;
    
end

model = Battery(paramobj);

switch schedulecase

  case 1

    % Schedule with two phases : activation and operation
    % 
    % Activation phase with exponentially increasing time step
    n = 25; 
    dt = []; 
    dt = [dt; repmat(0.5e-4, n, 1).*1.5.^[1:n]']; 
    % Operation phase with constant time step
    n = 40; 
    dt = [dt; repmat(2e-1*hour, n, 1)]; 
    
    % Time scaling can be adding using variable tfac
    times = [0; cumsum(dt)]*tfac; 
    
  case 2

    % Schedule used in activation test 
    n = 10;
    dt = rampupTimesteps(1.5*tup, tup/n, 10);
    times = [0; cumsum(dt)]; 

  case 3
    
    % Schedule adjusted for 1D case
    dt1 = rampupTimesteps(0.1, 0.1, 5);
    dt2 = 0.1*hour*ones(30, 1);
    dt = [dt1; dt2];
    times = [0; cumsum(dt)]; 

  case 4

    % Schedule with two phases : activation and operation
    % 
    % Activation phase with exponentially increasing time step
    n = 5; 
    dt = []; 
    dt = [dt; repmat(0.5e-4, n, 1).*1.5.^[1:n]']; 
    % Operation phase with constant time step
    %n = 24; 
    %dt = [dt; dt(end).*2.^[1:n]']; 
    %dt = [dt; repmat(dt(end)*1.5, floor(n*1.5), 1)]; 
    
    % Time scaling can be adding using variable tfac
    times = [0; cumsum(dt)]*tfac; 

  case 5

    % Schedule with two phases : activation and operation
    % 
    % Activation phase with exponentially increasing time step
    n = 25; 
    dt = []; 
    dt = [dt; repmat(0.5e-4, n, 1).*1.5.^[1:n]']; 
    % Operation phase with constant time step
    n = 40; 
    dt = [dt; repmat(4e-1*hour, n, 1)]; 
    
    % Time scaling can be adding using variable tfac
    times = [0; cumsum(dt)]; 
    
end

tt = times(2 : end); 
initstate = model.setupInitialState(); 

step = struct('val', diff(times), 'control', ones(numel(tt), 1)); 

pe = 'PositiveElectrode';
cc = 'CurrentCollector';
stopFunc = @(model, state, state_prev) (state.(pe).(cc).E < 2.0); 

srcfunc = @(time) CurrentSource(time, tup, times(end), model.J); 

control = repmat(struct('src', srcfunc, 'stopFunction', stopFunc), 1, 1); 
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
time = cellfun(@(x) x.time, states); 

%% plot

doplotJ = true;
if doplotJ
    figure
    for i = 1 : numel(time)
        src(i) = srcfunc(time(i));
    end
    plot((time/hour), src, '*-')
    xlabel('time (hours)')
    title('J (sine rampup)')
end

figure
plot((time/hour), Enew, '*-')
title('Potential (E)')
xlabel('time (hours)')

%% more plotting (case dependent)

switch modelcase
  case '1D'
    % plot1D;
  case '2D'
    plotThermal(model,states);
  case '3D'
    plot3D;
    plot3Dconc;
    plot3Dphi;
end