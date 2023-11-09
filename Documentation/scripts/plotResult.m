function plotResult(output)

    states = output.states;

    for ind = 1 : numel(states)
        state = states{ind};
        
        time(ind) = states{ind}.time;
        E(ind)    = states{ind}.Control.E;
        
    end

    plot(time/hour, E);

    xlabel('time / h');
    ylabel('voltage / V');
    
end
