jsonstruct = parseBattmoJson('simulation.json');

output = runBattery(jsonstruct);

%% plotting

t = output.time;
E = output.E;

figure
plot(t/hour, E)
xlabel('Time / h')
xlabel('Voltage / V')


