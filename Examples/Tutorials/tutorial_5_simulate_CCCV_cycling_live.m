%% Tutorial 5 - Simulate CC-CV Cycling
%% Introduction
% In this tutorial, we will use a P2D model to simulate CC-CV cycling. After 
% completing this tutorial, you should have a working knowledge of:
%% 
% * How to define and modify cycling protocols in BattMo
%% 
% We'll use the same model from Tutorial 1.

jsonstruct = parseBattmoJson('Examples/jsondatafiles/sample_input.json');
%% 
% Parameters are defined in the JSON parameter file and parsed into the MATLAB 
% structure. Once the JSON file has been read into MATLAB as a jsonstruct, its 
% properties can be modified programmatically.
%% Explore the Control Definition
% Let's begin by reviewing the control protocol in BattMo, with the command:

disp(jsonstruct.Control)
%% 
% We see that the default control protocol is set to a constant current (galvanostatic) 
% discharge. To change to a CC-CV cycling protocol, we can use the command:

cccv_control_protocol = parseBattmoJson('cccv_control.json');
jsonstruct_modified = mergeStructs({cccv_control_protocol, jsonstruct});
%% 
% Now we can explore the modified control protocol definition with the command:

disp(jsonstruct_modified.Control)
%% 
% Let's run the simulation and plot the cell voltage curve.

% run the simulation
output = runBatteryJson(jsonstruct_modified);
%% 
% get the states

states = output.states;

% extract the time and voltage quantities
time = cellfun(@(state) state.time, states);
voltage = cellfun(@(state) state.('Control').E, states);
current = cellfun(@(state) state.('Control').I, states);

% calculate the capacity
capacity = time .* current;

% plot the discharge curve in the figure
plot(time/hour, voltage, '-', 'linewidth', 3)

% add plot annotations
xlabel('Capacity  /  mA \cdot h')
ylabel('Cell Voltage  /  V')
%% Summary
% In this tutorial, we explored how to modify material parameters in BattMo.