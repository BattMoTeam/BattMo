clear all
close all

% setup mrst modules
mrstModule add ad-core multimodel mrst-gui battery agmg mpfa

mrstVerbose off

% Value used in rampup function, see currentSource.
tup = 0.1;

paramobj = LithiumBatteryInputParams();

% Setup battery

modelcase = '1D';

switch modelcase

  case '1D'

    gen = BatteryGenerator1D();
    paramobj = gen.updateBatteryInputParams(paramobj);
    schedulecase = 3; 

  case '2D'

    gen = BatteryGenerator2D();
    paramobj = gen.updateBatteryInputParams(paramobj);
    paramobj.elyte.thermalConductivity = 1;
    paramobj.elyte.heatCapacity = 1;
    schedulecase = 1;

    tfac = 1; % used in schedule setup
  
  case '3D'

    gen = BatteryGenerator3D();
    fac=1;
    gen.facx=fac;
    gen.facy=fac;
    gen.facz=fac;
    gen = gen.applyResolutionFactors();
    schedulecase = 1;
    %schedulecase = 4;% for testing
    tfac = 40; % used in schedule setup
    
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
    n = 24; 
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
    dt1 = rampupTimesteps(0.1, 0.1, 5);
    dt2 = 3e3*ones(30, 1);
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
    nls.LinearSolver = LinearSolverBattery('method', 'agmg', 'verbosity', 1);
    nls.LinearSolver.tol = 1e-3;
    nls.verbose = 10
end
model.nonlinearTolerance = 1e-5; 
model.verbose = false;

% Run simulation

doprofiling = true;
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


%% 1D plot

if strcmp(modelcase, '1D')

    h = figure;
    
    ffields = {'phi', 'c', 'j', 'LiFlux'};

    mnames = {{'Electrolyte'}, ...
              {'PositiveElectrode','ElectrodeActiveComponent'}, ...
              {'NegativeElectrode','ElectrodeActiveComponent'}, ...
              {'NegativeElectrode','CurrentCollector'}, ...
              {'PositiveElectrode','CurrentCollector'}};    
    
    tM = max(states{end}.T);
    tm = min(states{1}.T);
    tM = tM + 1e-1*(tM - tm);
    tm = tm - 1e-1*(tM - tm);
    
    for i = 1 : numel(states)
        
        state = states{i}; 

        for k = 1 : numel(ffields)

            ffield = ffields{k};
            
            figure(h)
            subplot(2, 3, k)
            cla, hold on
            
            if (strcmp(ffield, 'phi') || strcmp(ffield, 'j'))
                inds = 1 : 5; 
            else
                inds = 1 : 3; 
            end
            
            for ind = inds
                
                mname = mnames{ind}; 
                submodel = model.getSubmodel(mname); 
                substate = model.getProp(state, mname); 
                
                if strcmp(ffield, 'c') && (ind == 1)
                    var = substate.cs{1}; 
                else
                    var = substate.(ffield); 
                end
                
                if (k < 3)
                    plot(submodel.G.cells.centroids(:, 1), var, '* - ')
                else
                    iface = all(submodel.G.faces.neighbors > 0, 2); 
                    plot(submodel.G.faces.centroids(iface, 1), var, '* - ')
                end
                
                subtitle(ffield)
                
            end
        end
        
        % plot temperature
        subplot(2, 3, 5);
        c = model.G.cells.centroids(:, 1);
        plot(model.G.cells.centroids(:, 1), state.T, '* - ');
        axis([min(c), max(c), tm, tM]);
        drawnow;
        pause(0.01);
        
    end
    
end