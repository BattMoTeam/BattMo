%% Equivalent Circuit Model Simulation

% Load data set
jsonstruct = jsondecode(fileread('exampleECM.json'));
inputparams = EquivalentCircuitModelInputParams(jsonstruct);

% Setup model
model = EquivalentCircuitModel(inputparams);

% Run Simulation
[t, U, I, SOC] = model.solve();


% Plotting
figure; 
tiledlayout(3, 1);

nexttile
plot(t, U, 'LineWidth', 2)
title('Tension de la Batterie')
xlabel('Temps (s)')
ylabel('Tension (V)')
grid on

nexttile
plot(t, I, 'LineWidth', 2, 'Color', 'r')
title('Courant Appliqué')
xlabel('Temps (s)')
ylabel('Courant (A)')
grid on

nexttile
plot(t, SOC, 'LineWidth', 2, 'Color', 'g')
title('SOC')
xlabel('Temps (s)')
ylabel('SOC')
grid on