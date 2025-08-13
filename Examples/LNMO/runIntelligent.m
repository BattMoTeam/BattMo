clear
close all
clc

jsonstruct = parseBattmoJson('params_initial.json');

% Control
jsonfilename = fullfile('Examples', 'JsonDataFiles', 'cccv_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);

jsonstruct_control.Control.lowerCutoffVoltage = 4;
jsonstruct_control.Control.upperCutoffVoltage = 4.8;
jsonstruct_control.Control.DRate              = 1;
jsonstruct_control.Control.CRate              = 1;
jsonstruct_control.Control.dIdtLimit          = 1E-4;
jsonstruct_control.Control.dEdtLimit          = 1E-2;
jsonstruct_control.Control.numberOfCycles     = 5;

jsonstruct = mergeJsonStructs({jsonstruct_control, jsonstruct});

jsonstruct.TimeStepping.numberOfTimeSteps = 1000;

% For this simulation, to resolve the control switc we need to increase the number of time step cut allowed.
jsonstruct.NonLinearSolver.maxTimestepCuts = 20;

output = runBatteryJson(jsonstruct, 'runSimulation', true);

% Define the vectors
I = output.I; % Vector containing current values
t = output.time; % Vector containing corresponding time values

% Initialize the capacity vector
capacity = zeros(size(I));

% Integrate the current over time to calculate capacity
for i = 2:length(I)
    % Accumulate the integral based on the sign of the current
    capacity(i) = capacity(i-1) + (I(i) + I(i-1)) / 2 * (t(i) - t(i-1));
end

% Ensure the capacity starts at 0 and reverses back to 0
% capacity = capacity - min(capacity); % Shift to start from 0
% capacity = max(capacity) - capacity; % Reverse back to 0 when fully charged


figure('Units', 'centimeters', 'Position', [0, 0, 29.7, 21]);
plot(capacity ./ 3600, output.E)
xlabel('Capacity  /  A \cdot h')
ylabel('Cell Voltage  /  V')
set(gca, 'FontSize', 18);

%plotContours()
% 
% states = output.states;
% time = output.time;
% x = output.model.grid.cells.centroids;
% 
% 
% % instantiate matrices
% numStates = numel(states);
% maxSize = max(cellfun(@(x) numel(x.Electrolyte.c), output.states));
% c_elyte = NaN(numStates, maxSize);
% 
% % Loop through each element in 'states' and extract 'c' values
% for i = 1:numStates
%     c_elyte(i, 1:numel(output.states{i}.Electrolyte.c)) = output.states{i}.Electrolyte.c;
% end
% 
% % plot the concentration values over space (x axis) and time (y axis)
% h = figure();
% contourf(x, time/hour, c_elyte, 20, 'LineWidth',0.1)
% cm = cmocean('curl', 'pivot', 1000);
% colormap(cm)
% xlabel('Position  /  µm')
% ylabel('Time  /  h')
% colorbar()
