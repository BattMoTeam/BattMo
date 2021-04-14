clear all
close all

% setup mrst modules
mrstModule add ad-core multimodel mrst-gui battery matlab_bgl
mrstVerbose off

% Value used in rampup function, see currentSource.
tup = 0.1;

paramobj = LithiumBatteryInputParams();

% setup battery

modelcase = '1D';

switch modelcase
  case '1D'
    gen = BatteryGenerator1D();
    schedulecase = 3; 
  case '2D'
    gen = BatteryGenerator2D();
    schedulecase = 1;
    tfac = 10; % used in schedule setup
  case '3D'
    gen = BatteryGenerator3D();
    schedulecase = 1;
    tfac = 40; % used in schedule setup
end
if(true)
fac=1/5
gen.sepnx=10*fac;
gen.nenx=10*fac;
gen.penx=10*fac;
gen.ccnenx=10*fac;
gen.ccpenx=10*fac;
end
paramobj = gen.updateBatteryInputParams(paramobj);
if(false)
fac_e=10*1e-6;
fac_p=10*1e-8;
paramobj.ne.eac.am.electronicConductivity=100*fac_e;
paramobj.ne.cc.EffectiveElectronicConductivity=100*fac_e;
paramobj.pe.eac.am.electronicConductivity=100*fac_e;
paramobj.pe.cc.EffectiveElectronicConductivity=100*fac_e;
paramobj.pe.eac.am.electronicConductivity=100*fac_e;
paramobj.elyte.conductivityFactor=fac_p;
end
model = Battery(paramobj);

switch schedulecase

  case 1

    % Schedule with two phases : activation and operation
    % 
    % Activation phase with exponentially increasing time step
    n = 4;%25; 
    dt = []; 
    dt = [dt; repmat(1e-4, n, 1).*1.5.^[1:n]']; 
    % Operation phase with constant time step
    n = 24; 
    %dt = [dt; dt(end).*2.^[1:n]']; 
    %dt = [dt; repmat(dt(end)*1.5, floor(n*1.5), 1)]; 
    
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
    dt2 = 10*5e3*ones(40, 1);
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

%nls.LinearSolver=LinearSolverBattery('method','iterative');
%nls.LinearSolver=LinearSolverBattery('method','direct');
%nls.LinearSolver=LinearSolverBattery('method','agmg');
% Run simulation
model.AutoDiffBackend=DiagonalAutoDiffBackend();
%%
%profile -detail builtin
profile off
profile on
[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule,...
                                                'OutputMinisteps', true,...
                                                'NonLinearSolver', nls); 
profile off
profile report
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
return

if(strcmp(modelcase,'1D'))
    %% 1D plot
    sind=[12,34]
    figure(1),clf,hold on
    
    ffields = {'phi','c','j','LiFlux'}
    
    for k=1:numel(ffields)
        mnames={{'Electrolyte'},{'PositiveElectrode','ElectrodeActiveComponent'},{'NegativeElectrode','ElectrodeActiveComponent'},...
            {'NegativeElectrode','CurrentCollector'},{'PositiveElectrode','CurrentCollector'}}
        subplot(2,2,k),cla,hold on
        if(strcmp(ffields{k},'phi') || strcmp(ffields{k},'j') )
            ind=1:5;
        else
            ind=1:3;
        end
        
        for kk=sind
            state = states{kk};
            for i=ind%
                %for i=[2,3,1]
                mname=mnames{i};
                submodel = model.getSubmodel(mname)
                substate = model.getProp(state,mname);
                %fname=snames{i};
                if(not(isempty(state)))
                    if(strcmp(ffields{k},'c') && i==1)
                        var=substate.cs;
                    else
                        var = substate.(ffields{k});
                    end
                    if(iscell(var))
                        var=var{1};
                    end
                    if(k<3)
                        plot(submodel.G.cells.centroids(:,1), var,'*-')
                        hold on
                    else
                        iface = all(submodel.G.faces.neighbors>0,2);
                        plot(submodel.G.faces.centroids(iface,1), var,'*-')
                        hold on 
                    end
                    subtitle(ffields{k})
                    %pause()
                end
            end
        end
    end
    %%
end
