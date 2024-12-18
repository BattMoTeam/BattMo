%% Gas Supply Layer
%
% Single gas supply layer where a mixture of two gases (here H2O and O2) is injected at the bottom and flows to the top
% (See illustration below).
%

%% json input data
%
% We load the input data that is given by json structure
%

filename = fullfile(battmoDir(), 'ProtonicMembrane', 'jsonfiles', 'gas_supply.json');
jsonstruct_material = parseBattmoJson(filename);

filename = fullfile(battmoDir(), 'ProtonicMembrane', 'jsonfiles', '2d-gas-layer-geometry.json');
jsonstruct_geometry = parseBattmoJson(filename);

jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                              jsonstruct_geometry});

%% Input parameter setup
%
% We setup the input parameter structure which will we be used to instantiate the model
%

inputparams = ProtonicMembraneGasSupplyInputParams(jsonstruct);

%%
% We setup the grid, which is done by calling the function :battmo:`setupProtonicMembraneGasLayerGrid`
[inputparams, gridGenerator] = setupProtonicMembraneGasLayerGrid(inputparams, jsonstruct);

%% Model setup
%
% We instantiate the model for the gas supply layer
%

model = ProtonicMembraneGasSupply(inputparams);

%%
% The model is equipped for simulation using the following command (this step may become unnecessary in future versions)

model = model.setupForSimulation();

%% Model Plot
%
% We plot the simulation grid
%
G = model.grid;

plotGrid(G);

%%
% The boundary conditions have been given by the input json file, see :battmofile:`gas_supply.json <ProtonicMembrane/jsonsfiles/gas_supply.json>`
%
coupTerm = model.couplingTerms{1};
plotFaces(G, coupTerm.couplingfaces, 'edgecolor', 'blue', 'linewidth', 3)

coupTerm = model.couplingTerms{2};
plotFaces(G, coupTerm.couplingfaces, 'edgecolor', 'red', 'linewidth', 3)


return

%% Setupinitial state
%
% The model provides us with a default initialization

initstate = model.setupInitialState();

%%
% We modify the initial state by changin the initial mass fractions

mfInit = 0.2;
initstate.massfractions{1}             = mfInit + 0*initstate.massfractions{1};
initstate.GasSupplyBc.massfractions{1} = mfInit + 0*initstate.GasSupplyBc.massfractions{1};

%% Schedule setup
%
% The schedule structure constains two fields :code:`step`, :code:`control`, which gives respectively the time steps and
% the control to be applied.

T  = 1e-2*second;
N  = 100;

dt = rampupTimesteps(T, T/N, 1);

step.val = dt;
step.control = ones(numel(step.val), 1);

control.src = []; % The control is taking care of by a dedicated submobel. Hence, the empty field here.

schedule = struct('control', control, 'step', step);

%% Simulation
%
% We start the simulation
[~, states, report] = simulateScheduleAD(initstate, model, schedule);


%%

figure
plotToolbar(model.grid, states);
% caxis([0.2, 0.4])
% uit = findobj(gcf, 'Tooltip', 'Freeze caxis');
% uit.State = 'on';



