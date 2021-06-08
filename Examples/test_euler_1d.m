clear all
close all

% setup mrst modules
mrstModule add ad-core mrst-gui battery
mrstVerbose off

modelcase = '2D';

switch modelcase
  case '1D'
    inputparams = BatteryInputParams1D();
    inputparams.I = 1;
    schedulecase = 3;
    cond=100*1e-9;
    inputparams.NegativeElectrode.ActiveMaterial.electricalConductivity=cond
    inputparams.PositiveElectrode.ActiveMaterial.electricalConductivity=cond
    inputparams.PositiveCurrentCollector.electricalConductivity=cond;
    inputparams.NegativeCurrentCollector.electricalConductivity=cond;
    %inputparams.Electrolyte.eps=
  case '2D'
    inputparams = BatteryInputParams2D();
    schedulecase = 1;
    %inputparams.I = 1e-4;
    tfac = 1; % used in schedule setup
  case '3D_1'
    inputparams = BatteryInputParams3D_1();
    schedulecase = 1;
    tfac = 1; % used in schedule setup
  case '3D_2'
    inputparams = BatteryInputParams3D_2();
    inputparams.I = 1e-4;
    schedulecase = 1;
    tfac = 40; % used in schedule setup
end

% setup Battery model using parameter inputs.
model = Battery(inputparams);

% Value used in rampup function, see currentSource.
tup = 0.1;

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

stopFunc = @(model, state, state_prev) (state.PositiveCurrentCollector.E < 2.0); 

srcfunc = @(time) CurrentSource(time, tup, times(end), model.I); 

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
nls.LinearSolver=LinearSolverBattery('method','iterative');
%nls.LinearSolver=LinearSolverBattery('method','direct');

% Run simulation
[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule,...
                                                'OutputMinisteps', true,...
                                                'NonLinearSolver', nls); 

%%  Process output

ind = cellfun(@(x) not(isempty(x)), states); 
Enew = cellfun(@(x) x.PositiveCurrentCollector.E, {states{ind}}); 
time = cellfun(@(x) x.time, {states{ind}}); 

%% plot

doplotI = true;
if doplotI
    figure
    for i = 1 : numel(time)
        src(i) = srcfunc(time(i));
    end
    plot((time/hour), src, '*-')
    xlabel('time (hours)')
    title('I (sine rampup)')
end

figure
plot((time/hour), Enew, '*-')
title('Potential (E)')
xlabel('time (hours)')

%% 1D plot
ind=11
figure(1),clf,hold on
subplot(3,2,1),cla,hold on
mnames={'Electrolyte','NegativeCurrentCollector','PositiveCurrentCollector','PositiveElectrode','NegativeElectrode'}
snames={{'Electrolyte','phi'},{'NegativeCurrentCollector','phi'},{'PositiveCurrentCollector','phi'},...
    {'PositiveElectrode','ActiveMaterial','phi'},{'NegativeElectrode','ActiveMaterial','phi'}};
%for k=1:5:numel(states)
for k=ind
    state = states{k};
    for i=1:numel(mnames)%
        %for i=[2,3,1]
        mname=mnames{i};
        fname=snames{i};
        if(not(isempty(state)))
            plot(model.(mname).G.cells.centroids(:,1), model.getProps(state,fname),'*-')
            hold on
        end
        %pause()
    end
end

%
%figure(2),cla,hold on
subplot(3,2,2),cla, hold on
mnames={'Electrolyte','PositiveElectrode','NegativeElectrode'}
snames={{'Electrolyte','cs'},{'PositiveElectrode','ActiveMaterial','Li'},{'NegativeElectrode','ActiveMaterial','Li'}};
%for k=1:5:numel(states)
for k=ind
    state = states{k};
    for i=1:numel(mnames)%
        %for i=[2,3,1]
        mname=mnames{i};
        fname=snames{i};
        if(not(isempty(state)))
            var = model.getProps(state,fname);
            if(iscell(var))
                var=var{1};
            end
            plot(model.(mname).G.cells.centroids(:,1), var,'*-')
            hold on
        end
        %pause()
    end
end


%
%figure(2),cla,hold on
subplot(3,2,3),cla, hold on
mnames={'Electrolyte','NegativeCurrentCollector','PositiveCurrentCollector','PositiveElectrode','NegativeElectrode'}
snames={{'Electrolyte','j'},{'NegativeCurrentCollector','j'},{'PositiveCurrentCollector','j'},...
    {'PositiveElectrode','j'},{'NegativeElectrode','j'}};
%for k=1:5:numel(states)
for k=ind
    state = states{k};
    for i=1:numel(mnames)%
        %for i=[2,3,1]
        mname=mnames{i};
        fname=snames{i};
        if(not(isempty(state)))
            var = model.getProps(state,fname);
            if(iscell(var))
                var=var{1};
            end
            iface = all(model.(mname).G.faces.neighbors>0,2);
            plot(model.(mname).G.faces.centroids(iface,1), var,'*-')
            hold on
        end
        %pause()
    end
end

%
%figure(2),cla,hold on
subplot(3,2,4),cla, hold on
mnames={'Electrolyte','PositiveElectrode','NegativeElectrode'}
snames={{'Electrolyte','LiFlux'},{'PositiveElectrode','LiFlux'},{'NegativeElectrode','LiFlux'}};
%for k=1:5:numel(states)
for k=ind
    state = states{k};
    for i=1:numel(mnames)%
        %for i=[2,3,1]
        mname=mnames{i};
        fname=snames{i};
        if(not(isempty(state)))
            var = model.getProps(state,fname);
            if(iscell(var))
                var=var{1};
            end
            iface = all(model.(mname).G.faces.neighbors>0,2);
            plot(model.(mname).G.faces.centroids(iface,1), var,'*-')
            hold on
        end
        %pause()
    end
end



%
%figure(2),cla,hold on
subplot(3,2,5),cla, hold on
mnames={'Electrolyte','PositiveElectrode','NegativeElectrode'}
snames={{'Electrolyte','eSource'},{'PositiveElectrode','eSource'},{'NegativeElectrode','eSource'}};
%for k=1:5:numel(states)
for k=ind
    state = states{k};
    for i=1:numel(mnames)%
        %for i=[2,3,1]
        mname=mnames{i};
        fname=snames{i};
        if(not(isempty(state)))
            var = model.getProps(state,fname);
            if(iscell(var))
                var=var{1};
            end
            plot(model.(mname).G.cells.centroids(:,1), var,'*-')
            hold on
        end
        %pause()
    end
end
%
%figure(2),cla,hold on
subplot(3,2,6),cla, hold on
mnames={'PositiveElectrode','NegativeElectrode'}
snames={{'PositiveElectrode','eta'},{'NegativeElectrode','eta'}};
%for k=1:5:numel(states)
for k=ind
    state = states{k};
    for i=1:numel(mnames)%
        %for i=[2,3,1]
        mname=mnames{i};
        fname=snames{i};
        if(not(isempty(state)))
            var = model.getProps(state,fname);
            if(iscell(var))
                var=var{1};
            end
            plot(model.(mname).G.cells.centroids(:,1), var,'*-')
            hold on
        end
        %pause()
    end
end