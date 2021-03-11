% _2D_test_case

clear all
close all

%% Add MRST module
mrstModule add ad-core multimodel mrst-gui
mrstVerbose on

model = BatteryModel3DSimpleextra();
% model = BatteryModel();
% model.verbose = true;
model.J = 0.1; 

%% plot of the computational graph

doplotgraph = false;
if doplotgraph
    g = setupGraph(model); 
    figure
    plot(g)
end

%% run simulation
% profile - detail builtin
doode15i = false; 

if doode15i
    profile off
    profile on
    [t, y] = model.p2d(); 
    profile off
    profile report

    initstate = icp2d(model); 
    model = setupFV(model, initstate); 
    sl = model.fv.slots; 
    clear state E
    for iy = 1 : size(y, 1)
        yy = y(iy, :)'; 
        states_elyte{iy}.Li = yy(sl{1}); 
        states_elyte{iy}.phi = yy(sl{2}); 
        states_ne{iy}.Li = yy(sl{3}); 
        states_ne{iy}.phi = yy(sl{4}); 
        states_pe{iy}.Li = yy(sl{5}); 
        states_pe{iy}.phi = yy(sl{6}); 
        states_ccne{iy}.phi = yy(sl{7}); 
        states_ccpe{iy}.phi = yy(sl{8}); 
        E(iy) = yy(sl{9}); 
    end

    plot(t/hour, E)
    return
end

%% run the same with euler

caseno = 2; 

switch caseno
    
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

%% 

% times = linspace(0, 1e-3, 100)'; 
tt = times(2:end); 
initstate = icp2d(model); 
step = struct('val', diff(times), 'control', ones(numel(tt), 1)); 

src = nan(numel(tt), 1); 
for i = 1:numel(tt)
    src(i) = currentSource(tt(i), 0.1, 86400, model.J); 
end

stopFunc = @(model, state, state_prev) (state.ccpe.E < 2.0); 
srcfunc = @(time) currentSource(time, 0.1, 86400, model.J); 

control = repmat(struct('src', srcfunc, 'stopFunction', stopFunc), 1, 1); 
schedule = struct('control', control, 'step', step); 
nls = NonLinearSolver(); 

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

initstate.wellSol = []; 
% profile - detail builtin
% afterStepFn = @(model, states, reports, solver, schedule, simtime)....
% deal(model, states, reports, solver, states{end}.ccpe.E<2); 

if(false)
    % nls.timeStepSelector = IterationCountTimeStepSelector('targetIterationCount', 5); 
    nls.timeStepSelector = StateChangeTimeStepSelector(); 
    nls.timeStepSelector.targetProps = {{'ccpe', 'E'}}; 
    nls.timeStepSelector.targetChangeAbs = [0.1]; 
    nls.timeStepSelector.targetChangeRel = [1e9]; 
end
%nls.timeStepSelector = SimpleTimeStepSelector()
nls.maxIterations = 10; 
nls.errorOnFailure = false; 
model.nonlinearTolerance = 1e-4;
[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule,...
                                                'OutputMinisteps', true,...
                                                'NonLinearSolver', nls); 

%% 
ind = cellfun(@(x) not(isempty(x)), states); 
Enew = cellfun(@(x) x.ccpe.E, {states{ind}}); 
time = cellfun(@(x) x.time, {states{ind}}); 

figure
plot(log(time/hour), Enew, '*')
%plot((time/hour), Enew, '*')
title('Potential (E)')
xlabel('time (hours)')
%figure
%plotGrid(model.G)
return

ss = nan(size(t)); 
for i = 1:numel(t)
    ss(i) = schedule.control.src(t(i)); 
end
%% plotting 

% We set up the structure fv in this "hacky" way for the moment (it will be cleaned up in the future)
initstate = icp2d(model); 
model = setupFV(model, initstate); 
sl = model.fv.slots; 

clear state E
for iy = 1 : size(y, 1)
    yy = y(iy, :)'; 
    states_elyte{iy}.Li = yy(sl{1}); 
    states_elyte{iy}.phi = yy(sl{2}); 
    states_ne{iy}.Li = yy(sl{3}); 
    states_ne{iy}.phi = yy(sl{4}); 
    states_pe{iy}.Li = yy(sl{5}); 
    states_pe{iy}.phi = yy(sl{6}); 
    states_ccne{iy}.phi = yy(sl{7}); 
    states_ccpe{iy}.phi = yy(sl{8}); 
    E(iy) = yy(sl{9}); 
end

ind = cellfun(@(x) not(isempty(x)), states); 
% dt = schedule.step.val(ind); 

Enew = cellfun(@(x) x.ccpe.E, {states{ind}}); 
time = cellfun(@(x) x.time, {states{ind}}); 
figure
plot(time/hour, Enew, '*', t/hour, E)
title('Potential (E)')
xlabel('time (hours)')
% return
%% 
figure
plot(log(time/hour), Enew, '*', log(t/hour), E, ' - ')
title('Potential (E)')
xlabel('time (hours)')




