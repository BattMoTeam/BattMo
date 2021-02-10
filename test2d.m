%_2D_test_case

clear all
close all

%% Add MRST module
mrstModule add ad-core multimodel

model = BatteryModel();
% adminmodel = AdminModel();
% model = model.setupAdminModel(adminmodel);
model.J = 0.1;

g = setupGraph(model);

figure
plot(g)

[t, y] = model.p2d();

%% 

% We set up the structure fv in this "hacky" way for the moment (it will be cleaned up in the future)
initstate = icp2d(model);
model = setupFV(model, initstate);
fv = model.fv;

for iy = 1 : size(y, 1)
    yy = y(iy, :)';
    E(iy) = yy(model.fv.slots{end});
end

figure
plot(t/hour, E)
title('Potential (E)')
xlabel('time (hours)')

return
%%
n=30
dt=repmat(1e-3,n,1).*1.5.^[1:n]'
dt=[dt;repmat(dt(end),n*1.3,1)];
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
%%
tic;
[wellSols, states, report]  = simulateScheduleAD(initstate, model, schedule)
toc;
%% plot of each component

compnames = obj.componentnames;

allstates = {};
Gs = {};

for icn = 1 : numel(compnames)
    
    compname = compnames{icn};
    comp = obj.(compname);

    Gs{icn} = comp.Grid;
    
    allstates{icn}= {};
    
    for iy = 1 : size(y, 1)
        yy = y(iy, :)';
        varname = 'phi';
        fullvarname = sprintf('%s_%s', compname, varname);
        state.(varname) = yy(fv.getSlot(fullvarname));
        if ismember(compname, {'ne', 'pe', 'elyte'})
            varname = 'Li';
            fullvarname = sprintf('%s_%s', compname, varname);
            state.(varname) = yy(fv.getSlot(fullvarname));
            if strcmp(compname, 'elyte')
               state.(varname) = state.(varname)./obj.elyte.eps;
            end
        end
        allstates{icn}{iy} = state;
    end
    
    figure(icn)
    plotToolbar(Gs{icn}, allstates{icn});
    title(compname);
    
end

%% Combined plot for the positive electrod and current collector

% compnames = {'ne', 'ccne', 'pe', 'ccpe'};
% compnames = {'ne', 'ccne'};
compnames = {'pe', 'ccpe'};

states = {};
G = obj.G;
nc = G.cells.num;

for iy = 1 : size(y, 1)

    phi = nan(nc, 1);
    yy = y(iy, :)';
    
    for icn = 1 : numel(compnames)
    
        compname = compnames{icn};
        comp = obj.(compname);
        varname = 'phi';
        fullvarname = sprintf('%s_%s', compname, varname);
        philoc = yy(fv.getSlot(fullvarname));
        
        cellmap = comp.Grid.mappings.cellmap;
        
        phi(cellmap) = philoc;
        
    end

    states{iy}.phi = phi;
    
end

%%

figure
plotToolbar(G, states);
title(join(compnames, ' and '));


