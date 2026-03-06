%% Simulation Output Structure
% 

% jsonstruct = parseBattmoJson('Examples/Documentation/jsonfiles/explore_output_example.json');
% output = runBatteryJson(jsonstruct);

%% Output structure overview
%
% The output structure contains the results of the simulation, the simulation model and the simulation setup. We have in particular the folling
% 
% * |output.model| : The simulation model
% * |output.simsetup| : The simulation setup
% * |output.jsonstruct| : The json input file
% * |output.time| : Time structure of the results
% * |output.E| : Voltage for each time step
% * |output.I| : Current for each time step
% * |output.states| : The states structure with the internal state of the battery for each time step
%

%% Current and Voltage
%
% The readily available results are the current and voltage with the corresponding time structure
%

time = output.time;
E    = output.E;
I    = output.I;

figure
subplot(2,1,1)
plot(time/hour, E);
grid on
xlabel 'time  / h';
ylabel 'potential  / V';

subplot(2,1,2)
plot(time/hour, I);
grid on
xlabel 'time  / h';
ylabel 'Current  / I';


%% Simulation Model
%
% The simulation model is an output of the simulation. 

model = output.model;

%% 
%
% We can use it to plot the battery model, using the dedicated function |plotBatteryGrid|.
%

plotBatteryGrid(model);

%% Battery States
%
% The |states| gives us insight into the internal state of the battery. 
%

states = output.states;

%%
% They are stored in a cell array with the same time structure as the current and voltage. They follow the same
% architecture as the battery model, see <https://battmo.org/BattMo/architecture.html here>
%
%
% Let us consider the state computed at the last time step. 

state = states{end}

%%
%
% We obtain the following variables. The variables (except time and control type) are the _state_ or _primary_
% variables, which correspond to the unknown in the equations we solve. It is possible to recover also the intermediate
% variables as illustrated below.
% 
% * |state.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface|
% * |state.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion.c|
% * |state.NegativeElectrode.Coating.phi|
% * |state.NegativeElectrode.CurrentCollector.phi|
% * |state.PositiveElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface|
% * |state.PositiveElectrode.Coating.ActiveMaterial.SolidDiffusion.c|
% * |state.PositiveElectrode.Coating.phi|
% * |state.PositiveElectrode.CurrentCollector.scaledDeltaPhi|
% * |state.Electrolyte.phi|
% * |state.Electrolyte.c|
% * |state.ThermalModel.T|
% * |state.Control.E|
% * |state.Control.I|
% * |state.Control.ctrlType|
% * |state.time|
%

%%
% Let us plot the concentration in the electrolyte. The units are always defaulted to SI units. Note also the difference
% of scale for the different spatial directions.
%

figure
plotCellData(model.Electrolyte.grid, state.Electrolyte.c);
colorbar
view([50, 40]);

%% Simulation Setup
%
% The simulation setup contains all the information to re-run a simulation. It means the model, of course, but also the initial state, the schedule (time stepping choices) and the solver parameters.
%

simsetup = output.simsetup();
states = simsetup.run()

%% Json Input
%
% The |jsonstruct| that is returned in the output is the same as the one we sent as input but it is augmented with all the **default values**
%
%
