clear all
close all

% setup mrst modules
mrstModule add ad-core multimodel mrst-gui battery
mrstVerbose off

modelcase = '3D_2';

switch modelcase
  case '1D'
    inputparams = BatteryInputParams1D();
    inputparams.J = 1;
    schedulecase = 3; 
  case '2D'
    inputparams = BatteryInputParams2D();
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

% setup Battery model using parameter inputs.
model = BatteryModelSimple(inputparams);

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
    times = [0; cumsum(dt)]*tfac; 
    
    times(end); 
    
  case 2

    % Schedule with constant time steping
    n = 100;
    dt = repmat(1e-3, n, 1);
    times = [0; cumsum(dt)]; 
    times(end) 

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

stopFunc = @(model, state, state_prev) (state.ccpe.E < 2.0); 

srcfunc = @(time) currentSource(time, 0.1, times(end), model.J); 

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
Enew = cellfun(@(x) x.ccpe.E, {states{ind}}); 
time = cellfun(@(x) x.time, {states{ind}}); 

figure
%plot(log((time/hour)), Enew, '*')
plot((time/hour), Enew, '*')
title('Potential (E)')
xlabel('time (hours)')
