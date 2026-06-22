%% Equivalent Circuit Model Simulation

% Load data set

params = createParametersECM();
inputparams = EquivalentCircuitModelInputParams(params);

% Setup model
model = EquivalentCircuitModel(inputparams);

% Run Simulation
[t, U, I, SOC] = model.solve();


% Plotting
figure(Name='Highlighting the difference in resistance depending on the state of charge'); 
tiledlayout(3, 1);

nexttile
plot(t, U, 'LineWidth', 2)
title('Battery Voltage')
xlabel('Time /s')
ylabel('Voltage /V')
grid on

nexttile
plot(t, I, 'LineWidth', 2, 'Color', 'r')
title('Applied Current')
xlabel('Time /s')
ylabel('Current /A')
grid on

nexttile
plot(t, SOC, 'LineWidth', 2, 'Color', 'g')
title('SOC')
xlabel('Time /s')
ylabel('SOC')
grid on
