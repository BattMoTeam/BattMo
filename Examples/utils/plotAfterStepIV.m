function  [model, states, reports, solver, ok] = plotAfterStepIV(model, states, reports, solver, schedule, simtime)
    
    %% Process output and recover the output voltage and current from the output states.
ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
Enew = cellfun(@(x) x.Control.E, states); 
Inew = cellfun(@(x) x.Control.I, states);
Tmax = cellfun(@(x) max(x.ThermalModel.T), states);
[SOCN, SOCP] =  cellfun(@(x) model.calculateSOC(x), states);
time = cellfun(@(x) x.time, states);
dt = diff([0,time]);
SOC = cumsum(Inew.*dt)
figure(33),subplot(2,2,1),hold on
plot(time,Inew,'*-')
subplot(2,2,2),hold on
plot(SOC,Enew,'*-')
subplot(2,2,3),hold on
plot(time,Enew.*Inew.*dt,'*-')
subplot(2,2,4),hold on
SOC = cumsum(Inew.*dt)
plot(SOC,cumsum(Enew.*Inew.*dt),'*-')
ok = true;
figure(44),hold on
subplot(2,1,1)
plot(time,Tmax,'*-')
subplot(2,1,1)
plot(time,Tmax,'*-')
subplot(2,1,2)
plot(time,cumsum(Inew.*dt),'*-')
end