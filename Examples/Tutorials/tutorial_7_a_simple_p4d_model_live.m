%% Tutorial 7 - A Simple P4D Simulation
%% Introduction
% In this tutorial, we will use a simple P4D simulation to explore the effects 
% of battery cell architecture. After completing this tutorial, you should have 
% a working knowledge of:
%% 
% * How to setup and run a P4D simulation of a single cell in BattMo
% * Advanced usage of combining model descriptions from multiple sources
%% Construct the Model from Different Sources
% Let's say that you parameterize some cell materials and put that data in a 
% JSON file. Then you get some description of your cell geometry, and you put 
% that data in another JSON file. Now you want to combine those descriptions into 
% a single model and simulate it using some pre-defined contol protocol and simulation 
% settings in other files. How can we combine those easily, without having to 
% do any recoding?
% 
% To do that, we can simply make use of the <https://github.com/BattMoTeam/BattMo/blob/main/Utilities/JsonUtils/mergeStructs.m 
% mergeStructs> function in *BattMo* (an other examle is presented  <https://battmoteam.github.io/BattMo/mergejsonstruct.html#merging-parameters 
% here>.
% 
% First, let’s define our cell materials. We have provided a JSON file that 
% contains material properties for a NMC and Graphite active materials, which 
% we can parse as a *BattMo* structure:

% parse material definitions as a BattMo structure
jsonfilename = 'ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json';
jsonstruct_material = parseBattmoJson(jsonfilename);
%% 
% Next, we have defined the cell geometry properties in a separate JSON file 
% that we can also parse into *BattMo*:

% parse cell geometry specifications as a BattMo structure
jsonfilename = 'Examples/jsondatafiles/geometry3d.json';
jsonstruct_geometry = parseBattmoJson(jsonfilename);
%% 
% Let's have a closer look at the cell geometry specification.

disp(jsonstruct_geometry.Geometry)
%% 
% Here we can see that the width and height dimensions of the cell are defined, 
% along with the number of discretizations in each direction, and a case description. 
% The case sets the type of simulation to be performed. Here it is set to '3D-demo', 
% to BattMo knows to setp a P4D mesh.
% 
% We can take the same approach for the remaining parameters, as shown below, 
% for the control protocol, simulation settings, and output settings:

% control protocol
jsonfilename = fullfile('Examples', 'jsondatafiles', 'cc_discharge_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);

% simulation settings
jsonfilename = fullfile('Examples', 'jsondatafiles', 'simulation_parameters.json');
jsonstruct_simparams = parseBattmoJson(jsonfilename);

% output settings
jsonfilename = fullfile('Examples', 'jsondatafiles', 'extra_output.json');
jsonstruct_output = parseBattmoJson(jsonfilename);
%% 
% Now, we can merge these parameter definitions into a single parameter set 
% and run the simulation:

% combine the parameter structures from the different sources into a single
% BattMo structure
jsonstruct = mergeStructs({jsonstruct_geometry , ...
    jsonstruct_material , ...
    jsonstruct_control  , ...
    jsonstruct_simparams, ...
    jsonstruct_output   , ...
    });
%% 
% we store the output as a cell array so we can compare results across different 
% simulation runs in this tutorial; here we instantiate an empty cell array

output = cell(2,1);

% run the simulation
output{1} = runBattery(jsonstruct);
%% Visualize the Results
% We plot the model using <https://github.com/BattMoTeam/BattMo/blob/main/Utilities/Visualization/plotBatteryGrid.m 
% plotBatteryGrid> (note that the different axis are scaled differently)

% create a shorthand variable for the model
model = output{1}.model;

% use the plotBatteryGrid function to show the grid
plotBatteryGrid(model)

% make the axis tight and set the camera viewing angle
axis tight
view(45,45)
%% 
% We find a extensive set of plotting functions in <https://www.sintef.no/Projectweb/MRST/ 
% MRST>. You may be interested to have a look at the <https://www.sintef.no/projectweb/mrst/documentation/tutorials/visualization-tutorial/ 
% Visualization Tutorial>. Let us use the <https://github.com/SINTEF-AppliedCompSci/MRST/blob/main/core/plotting/plotGrid.m 
% plotGrid> and <https://github.com/SINTEF-AppliedCompSci/MRST/blob/main/core/plotting/plotCellData.m 
% plotCellData> to plot the surface particle concentrations in both electrode 
% at a given time step.

% set the timestep we want to visualize
timestep = 20;

% get the state of the simulation at the given timestep
state = output{1}.states{timestep};

% create a new figure
figure()

% plot the surface concentration of lithium in the negative electrode active material
plotCellData(model.NegativeElectrode.Coating.grid, state.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface/(mol/litre))

% plot the surface concentration of lithium in the positive electrode active material
plotCellData(model.PositiveElectrode.Coating.grid, state.PositiveElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface/(mol/litre))

% add a colorbar
colorbar()

% make the axis tight and set the camera viewing angle
axis tight
view(45,45)

% add plot annotations
title('Active Material Surface Lithium Concentration  /  mol \cdot L^{-1}');
%% Compare with a P2D Simulation

% change the setup of the model to consider a P2D case
jsonstruct.Geometry.case = '1D';
jsonstruct.Geometry.faceArea = jsonstruct.Geometry.width * jsonstruct.Geometry.height;

% change the rate
jsonstruct.Control.CRate = 1;

% update the total time of the simulation
jsonstruct.TimeStepping.totalTime = (1./jsonstruct.Control.CRate) .* 3600 .* 1.1;

% run the simulation
output{2} = runBattery(jsonstruct);
%% 
% In this case, MATLAB sends many warnings about the ill-conditionness of the 
% system. The ill-conditionness appears to come mainly from the very high ratio 
% between the electronic conductivity of the current collectors and the other 
% components.

% get the states from the P2D model
states_P2D = output{2}.states;

% extract the time and voltage quantities
time_P2D = cellfun(@(state) state.time, states_P2D);
voltage_P2D = cellfun(@(state) state.('Control').E, states_P2D);

% get the states from the P4D model
states_P4D = output{1}.states;

% extract the time and voltage quantities
time_P4D = cellfun(@(state) state.time, states_P4D);
voltage_P4D = cellfun(@(state) state.('Control').E, states_P4D);

% plot the discharge curves together in a new figure
figure();
plot((time_P2D/hour), voltage_P2D, '-', 'linewidth', 3)
hold on
plot((time_P4D/hour), voltage_P4D, '-', 'linewidth', 3)
xlabel('Time  /  h')
ylabel('Cell Voltage  /  V')
legend('P2D', 'P4D')
%% Summary
% In this tutorial, we learned how to create a simple P4D simulation in BattMo. 
% First, we explored how to combine parameter sets coming from a handful of different 
% files into a single coherent BattMo model description. Then we had a closer 
% look into the Geometry description to see how BattMo knows how to setup a P4D 
% or P2D type model. After running the P4D simulation, we learned how to visualize 
% simulation results on a 3D grid. For comparison, we ran the same model in a 
% P2D configuration and plotted the discharge curves together. This showed that 
% the results can diverge somewhat due to the effects of the tabs and non-ideal 
% transport in the electrode plane. These results show that P4D models can yield 
% important insight that may be lost in the averaged approach of P2D models.