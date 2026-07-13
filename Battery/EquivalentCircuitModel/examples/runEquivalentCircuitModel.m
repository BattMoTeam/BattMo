jsonstruct = createParametersECM();
inputparams =  EquivalentCircuitModelInputParams(jsonstruct);

model =  EquivalentCircuitModel(inputparams);

close all
cgit = model.cgit;
cgit.plot;


[t, U, I, SOC] = model.solve();

%%

figure
tiledlayout(1, 3)
nexttile
plot(t/hour, I)
title('Current')
xlabel('Time / hour')
ylabel('Current / A');
nexttile
plot(t/hour, U)
title('Voltage')
xlabel('Time / hour')
ylabel('Voltage / V');
nexttile
plot(t/hour, SOC)
title('SOC')
xlabel('Time / hour')
ylabel('SOC / -');
