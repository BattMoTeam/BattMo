clear all
close all

% setup mrst modules
mrstModule add ad-core multimodel mrst-gui battery
mrstVerbose off

modelcase = '2D';

switch modelcase
  case '2D'
    inputparams = BatteryInputParams2D();
  case '3D'
    inputparams = BatteryInputParams3D();
end

% setup Battery model using parameter inputs.
model = BatteryModelSimple(inputparams);

schedulecase = 2; 

switch schedulecase
    
  case 1
    dt = []; 
    n = 10
    dt = [dt; repmat(1e-3, n, 1).*1.5.^[1:n]']; 
    % n = 13
    % dt = [dt; dt(end).*2.^[1:n]']
    dt = [dt; repmat(dt(end), floor(n*1.5), 1)]; 
    times = [0; cumsum(dt)]; 
    times(end)
    
  case 2

    n = 15; 
    dt = []; 
    dt = [dt; repmat(1e-4, n, 1).*1.5.^[1:n]']; 
    n = 13; 
    dt = [dt; dt(end).*2.^[1:n]']; 
    dt = [dt; repmat(dt(end)*1.5, floor(n*1.5), 1)]; 
    times = [0; cumsum(dt)]; 
    times(end); 
    
  case 3
    dt = []; 
    % n = 37;
    % dt = [dt; repmat(1e-4, n, 1).*1.1.^[1:n]']; 
    n = 1;
    dt = [dt; repmat(1e-3, n, 1)];
    % n = 13
    % dt = [dt; dt(end).*2.^[1:n]']
    % dt = [dt; repmat(dt(end), floor(n*1.5), 1)]; 
    times = [0; cumsum(dt)]; 
    times(end) 
    
end

tt = times(2 : end); 
initstate = model.setupInitialState(); 
step = struct('val', diff(times), 'control', ones(numel(tt), 1)); 

src = nan(numel(tt), 1); 
for i = 1:numel(tt)
    src(i) = currentSource(tt(i), 0.1, 86400, model.J); 
end

stopFunc = @(model, state, state_prev) (state.ccpe.E < 2.0); 
srcfunc = @(time) currentSource(time, 0.1, 86400, model.J); 

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
plot((time/hour), Enew, '*')
title('Potential (E)')
xlabel('time (hours)')
