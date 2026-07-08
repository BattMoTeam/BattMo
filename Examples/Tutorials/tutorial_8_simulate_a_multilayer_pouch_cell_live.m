%% Tutorial 8 - Simulate a Multilayer Pouch Cell
%% Introduction
% In this tutorial, we simulate a multilayer pouch cell. We use the same material 
% property as in the other tutorials

jsonstruct_material = parseBattmoJson('Examples/jsondatafiles/sample_input.json');
%% 
% Next, we load and parse a json file where we have chosen some parameters for 
% the multilayer pouch domain. Note that all the parameters are described in a 
% json schema, see <https://github.com/BattMoTeam/BattMo/blob/main/Utilities/JsonSchemas/Geometry.schema.json 
% Geometry.schema.json>, even if the simplest way to proceed is to start with 
% an example, in this case given by <https://github.com/BattMoTeam/BattMo/blob/main/Examples/JsonDataFiles/geometryMultiLayerPouch.json 
% geometryMultiLayerPouch.json>.

jsonfilename = 'Examples/JsonDataFiles/geometryMultiLayerPouch.json';
jsonstruct_geometry = parseBattmoJson(jsonfilename);
%% 
% We use <https://github.com/BattMoTeam/BattMo/blob/main/Utilities/JsonUtils/FlatStructViewer.m 
% FlatStructViewer.m> to flatten the json structure and print it to screen. We 
% can see that, in this example, we use 5 layers and two different lengths for 
% the tabs (height value). At the moment, the two tabs share the same width. Implementing 
% a separate width for each tab would require to modify the grid generator for 
% this geometry. It is more a developper work but is definitely not out of reach.

fsv = flattenStruct(jsonstruct_geometry)
fsv.print();
%% 
% We load and parse the control protocol

jsonfilename = fullfile('Examples', 'jsondatafiles', 'cc_discharge_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);
%% 
% We load and parse the simulation settings. This is optional. Typically, reasonable 
% choices are made by default.

jsonfilename = fullfile('Examples', 'jsondatafiles', 'simulation_parameters.json');
jsonstruct_simparams = parseBattmoJson(jsonfilename);
%% 
% Now, we can merge these parameter definitions into a single parameter set 
% to obtain a jsonstruct that has all the input needed by the simulator.

jsonstruct = mergeStructs({jsonstruct_geometry , ...
    jsonstruct_material , ...
    jsonstruct_control  , ...
    jsonstruct_simparams}, 'warn', false);
%% Setup the model for inspection
% When we run the simulation using function <https://github.com/BattMoTeam/BattMo/blob/main/Examples/JsonInput/runBattery.m 
% runBattery.m>, the model is setup. In the case where we want to setup the 
% model for inspection, prior to simulation, we can use the function <https://github.com/BattMoTeam/BattMo/blob/main/Utilities/JsonUtils/setupModelFromJson.m 
% setupModelFromJson.m>

model = setupModelFromJson(jsonstruct);
%% 
% We use the <https://github.com/BattMoTeam/BattMo/blob/main/Utilities/Visualization/plotBatteryGrid.m 
% plotBatteryGrid.m> function to show the grid

plotBatteryGrid(model)
% make the axis tight and set the camera viewing angle
axis tight
view(45,45)
%% Run the simulation

output = runBattery(jsonstruct);
%% Visualize the Results
% extract the time and voltage quantities

states = output.states;

time    = cellfun(@(state) state.time, states);
voltage = cellfun(@(state) state.('Control').E, states);
%% 
% We plot the discharge curves together in a new figure

figure();
plot((time/hour), voltage, '-', 'linewidth', 3)
xlabel('Time  /  h')
ylabel('Cell Voltage  /  V')
title('Voltage');
%% 
% For a given time step, we plot the concentration on the grid.

% Set the timestep we want to visualize
timestep = 20;

% get the state of the simulation at the given timestep
state = states{timestep};

% create a new figure
figure()

% plot the surface concentration of lithium in the negative electrode active material
plotCellData(model.NegativeElectrode.Coating.grid, state.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface/(mol/litre))

% plot the surface concentration of lithium in the positive electrode active material
plotCellData(model.PositiveElectrode.Coating.grid, state.PositiveElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface/(mol/litre))

title('Active Material Surface Lithium Concentration  /  mol \cdot L^{-1}');
% add a colorbar
colorbar()
view(45,45)