function  [model, states, reports, solver, ok] = plotAfterStepIV(model, states, reports, solver, schedule, simtime)
%% Process output and recover the output voltage and current from the output states.

    ind = cellfun(@(x) not(isempty(x)), states); 
    states = states(ind);

    E    = cellfun(@(x) x.Control.E, states); 
    I    = cellfun(@(x) x.Control.I, states);
    Tmax = cellfun(@(x) max(x.ThermalModel.T), states);
    time = cellfun(@(x) x.time, states);
    
    dt = diff([0, time]);
    SOC = cumsum(I.*dt);

    %%
    figure(33)
    
    subplot(2, 2, 1)
    hold on
    plot(time, I, '*-')

    subplot(2, 2, 2), 
    hold on
    plot(SOC, E, '*-')

    subplot(2, 2, 3), 
    hold on
    plot(time, E.*I.*dt, '*-')

    subplot(2, 2, 4), 
    hold on
    plot(SOC, cumsum(E.*I.*dt), '*-')

    %%
    figure(44), 
    hold on
    
    subplot(2, 1, 1)
    plot(time, Tmax, '*-')
    
    subplot(2, 1, 1)
    plot(time, Tmax, '*-')
    
    subplot(2, 1, 2)
    plot(time, cumsum(I.*dt), '*-')

    ok = true;
end