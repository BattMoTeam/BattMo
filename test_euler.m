%_2D_test_case

clear all
close all

%% Add MRST module
mrstModule add ad-core multimodel mrst-gui

model = BatteryModel();
model.J = 0.1;

%% plot of the computational graph

g = setupGraph(model);
figure
plot(g)

%% run simulation
%profile -detail builtin
profile off
profile on
[t, y] = model.p2d();
profile off
profile report

%% run the same with euler
n=30
dt=repmat(1e-3,n,1).*1.5.^[1:n]'
dt=[dt;repmat(dt(end),floor(n*1.3),1)];
times=[0;cumsum(dt)];
%times=linspace(0,1e-3,100)';
tt=times(2:end);
initstate = icp2d(model)
step=struct('val',diff(times),'control',ones(numel(tt),1));
%
src=nan(numel(tt),1);
for i=1:numel(tt)
    src(i)=currentSource(tt(i), 0.1, 86400, model.J);
end
%
control=repmat(struct('src',[]),numel(src),1);
for i=1:numel(src)
    control(i).src=src(i);
end
step.control=[1:numel(src)]';
schedule=struct('control',control, 'step',step)
model.nonlinearTolerance=1e-3;
initstate.wellSol=[];
%profile -detail builtin
%afterStepFn =  @(model, states,  reports, solver, schedule, simtime)....
%      deal(model,states,reports,solver,states{end}.ccpe.E<2);
stopFunc =@(model,state) state.ccpe.E<2;
profile off
profile on
[wellSols, states, report]  = simulateScheduleAD(initstate, model, schedule,'OutputMinisteps',false,'stopFunction',stopFunc)
profile off
profile viewer


%% plotting 

% We set up the structure fv in this "hacky" way for the moment (it will be cleaned up in the future)
initstate = icp2d(model);
model = setupFV(model, initstate);
sl = model.fv.slots;

clear state E
for iy = 1 : size(y, 1)
    yy = y(iy, :)';
    states_elyte{iy}.Li  = yy(sl{1});
    states_elyte{iy}.phi = yy(sl{2});
    states_ne{iy}.Li     = yy(sl{3});
    states_ne{iy}.phi    = yy(sl{4});
    states_pe{iy}.Li     = yy(sl{5});
    states_pe{iy}.phi    = yy(sl{6});
    states_ccne{iy}.phi  = yy(sl{7});
    states_ccpe{iy}.phi  = yy(sl{8});
    E(iy) = yy(sl{9});
end

ind = cellfun(@(x) not(isempty(x)),states);
%dt = schedule.step.val(ind);

Enew= cellfun(@(x) x.ccpe.E,{states{ind}});
time = cellfun(@(x) x.t,{states{ind}});
figure
plot(time/hour, Enew,'*',t/hour, E)
title('Potential (E)')
xlabel('time (hours)')

%%
figure
plot(log(time/hour), Enew,'*',log(t/hour), E,'*')
title('Potential (E)')
xlabel('time (hours)')




