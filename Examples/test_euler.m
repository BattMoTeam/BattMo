clear all
close all

% setup mrst modules
mrstModule add ad-core multimodel mrst-gui battery matlab_bgl
mrstVerbose off

% Value used in rampup function, see currentSource.
tup = 0.1;

paramobj = LithiumBatteryInputParams();

% setup battery

modelcase = '2D';

switch modelcase
  case '1D'
    gen = BatteryGenerator1D();
    schedulecase = 3; 
  case '2D'
    gen = BatteryGenerator2D();
    schedulecase = 1;
    tfac = 1; % used in schedule setup
  case '3D_1'
    inputparams = BatteryInputParams3D_1();
    schedulecase = 1;
    tfac = 1; % used in schedule setup
  case '3D_2'
    inputparams = BatteryInputParams3D_2();
    inputparams.J = 1e-4;
    schedulecase = 1;
    tfac = 40; % used in schedule setup
end

paramobj = gen.updateBatteryInputParams(paramobj);

model = Battery(paramobj);

switch schedulecase

  case 1

    % Schedule with two phases : activation and operation
    % 
    % Activation phase with exponentially increasing time step
    n = 15; 
    dt = []; 
    dt = [dt; repmat(1e-4, n, 1).*1.5.^[1:n]']; 
    % Operation phase with constant time step
    n = 13; 
    dt = [dt; dt(end).*2.^[1:n]']; 
    dt = [dt; repmat(dt(end)*1.5, floor(n*1.5), 1)]; 
    
    % Time scaling can be adding using variable tfac
    times = [0; cumsum(dt)]; 
    
  case 2

    % Schedule used in activation test 
    n = 10;
    dt = rampupTimesteps(1.5*tup, tup/n, 10);
    times = [0; cumsum(dt)]; 

  case 3
    
    % Schedule adjusted for 1D case
    dt1 = rampupTimesteps(0.1, 0.1, 10);
    dt2 = 5e3*ones(40, 1);
    dt = [dt1; dt2];
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

% We can use a linear preconditioner (not yet fully tested, require extra setup)
useAMGCL = false;
if useAMGCL
    nls.LinearSolver = AMGCLSolverAD('verbose', true, 'write_params', true, 'preconditioner', 'relaxation', 'relaxation', ...
                                     'ilu0'); 
    nls.LinearSolver.verbose = true; 
    nls.LinearSolver.amgcl_setup.verbose
    nls.LinearSolver.amgcl_setup.verbose = true
    nls.LinearSolver = AGMGSolverAD(); 
    nls.LinearSolver.reduceToCell = true
    nls.LinearSolver.tolerance = 1e-5
end

% Run simulation
[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule,...
                                                'OutputMinisteps', true,...
                                                'NonLinearSolver', nls); 

%%  Process output

ind = cellfun(@(x) not(isempty(x)), states); 
Enew = cellfun(@(x) x.(pe).(cc).E, {states{ind}}); 
time = cellfun(@(x) x.time, {states{ind}}); 

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



