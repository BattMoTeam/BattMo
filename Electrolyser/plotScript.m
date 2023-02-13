time = cellfun(@(state) state.time, states);
E = cellfun(@(state) state.(oer).(ptl).E, states);

close all
figure
plot(time, E)


