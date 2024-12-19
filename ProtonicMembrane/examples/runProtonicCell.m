%% PEM electrolyser with Gas Supply
%

%% json input data
%
% We load the input data that is given by json structures. The physical properties of the supply gas layer is given in
% the json file :battmofile:`gas-supply-whole-cell.json <ProtonicMembrane/jsonfiles/gas-supply-whole-cell.json>`

clear jsonstruct_material
filename = 'ProtonicMembrane/jsonfiles/gas-supply-whole-cell.json';
jsonstruct_material.GasSupply = parseBattmoJson(filename);
filename = 'ProtonicMembrane/jsonfiles/protonicMembrane.json';
jsonstruct_material.Electrolyser = parseBattmoJson(filename);

jsonstruct_material = removeJsonStructFields(jsonstruct_material, ...
                                             {'Electrolyser', 'Electrolyte', 'Nx'}, ...
                                             {'Electrolyser', 'Electrolyte', 'xlength'});

filename = 'ProtonicMembrane/jsonfiles/2d-cell-geometry.json';
jsonstruct_geometry = parseBattmoJson(filename);

filename = 'ProtonicMembrane/jsonfiles/gas-supply-initialization.json';
jsonstruct_initialization.GasSupply = parseBattmoJson(filename);

jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                               jsonstruct_geometry, ...
                               jsonstruct_initialization});


%% Input parameter setup
%
% We setup the input parameter structure which will we be used to instantiate the model
%

inputparams = ProtonicMembraneCellInputParams(jsonstruct);

%%
%
% We setup the grid, which is done by calling the function :battmo:`setupProtonicMembraneGasLayerGrid`
%

[inputparams, gen] = setupProtonicMembraneCellGrid(inputparams, jsonstruct);

%% Model setup
%
%

model = ProtonicMembraneCell(inputparams);

%%
% The model is equipped for simulation using the following command (this step may become unnecessary in future versions)

model = model.setupForSimulation();

%% Grid plots
%

figure('position', [337, 757, 3068, 557])
hold on
plotGrid(model.grid)
plotGrid(model.Electrolyser.grid, 'facecolor', 'red')
plotGrid(model.GasSupply.grid, 'facecolor', 'blue')

plotGrid(model.GasSupply.grid, model.GasSupply.couplingTerms{1}.couplingcells);
plotGrid(model.GasSupply.grid, model.GasSupply.couplingTerms{2}.couplingcells);
plotGrid(model.Electrolyser.grid, model.Electrolyser.couplingTerms{1}.couplingcells(:, 2));
plotGrid(model.Electrolyser.grid, model.Electrolyser.couplingTerms{2}.couplingcells(:, 2));
plotGrid(model.GasSupply.grid, model.couplingTerm.couplingcells(:, 1) );

%% Setup initial state
%

initstate = model.setupInitialState(jsonstruct);

%% Setup schedule
%

schedule = model.setupSchedule(jsonstruct);

return

%% Setup nonlinear solver

ts = IterationCountTimeStepSelector('targetIterationCount', 5, ...
                                    'verbose', true);

nls = NonLinearSolver();
nls.timeStepSelector = ts;
nls.maxIterations    = 15;
nls.errorOnFailure   = false;
nls.verbose          = true;

model.nonlinearTolerance = 1e-5;
model.verbose = true;

%% Start simulation

dopack = true;
clearSimulation = true;

if dopack

    name = 'testwholecell';
    problem = packSimulationProblem(initstate, model, schedule, []             , ...
                                    'ExtraArguments', {'OutputMinisteps', true}, ...
                                    'Name'           , name                    , ...
                                    'NonLinearSolver', nls);
    if clearSimulation
        %% clear previously computed simulation
        clearPackedSimulatorOutput(problem, 'prompt', false);
    end
    simulatePackedProblem(problem);
    [globvars, states, reports] = getPackedSimulatorOutput(problem);
    
else
    
    [~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

end

%% plotting

close all

set(0, 'defaultlinelinewidth', 3);
set(0, 'defaultaxesfontsize', 15);

N = gen.nxElectrolyser;
xc = model.(elyser).(elyte).grid.cells.centroids(1 : N, 1);

state = states{end};

state = model.addVariables(state);

X = reshape(model.(elyser).(elyte).grid.cells.centroids(:, 1), N, [])/(milli*meter);
Y = reshape(model.(elyser).(elyte).grid.cells.centroids(:, 2), N, [])/(milli*meter);

figure
val = state.(elyser).(elyte).pi;
Z = reshape(val, N, []);
surf(X, Y, Z, 'edgecolor', 'none');
title('pi / V')
xlabel('x / mm')
ylabel('y / mm')
view(45, 31)
colorbar

figure
val = state.(elyser).(elyte).pi - state.(elyser).(elyte).phi;
Z = reshape(val, N, []);
surf(X, Y, Z, 'edgecolor', 'none');
title('E / V')
xlabel('x / mm')
ylabel('y / mm')
view(45, 31)
colorbar

figure
val = state.(elyser).(elyte).phi;
Z = reshape(val, N, []);
surf(X, Y, Z, 'edgecolor', 'none');
title('phi / V')
xlabel('x / mm')
ylabel('y / mm')
view(73, 12)
colorbar

N = gen.nxGasSupply;

X = reshape(model.(gs).grid.cells.centroids(:, 1), N, [])/(milli*meter);
Y = reshape(model.(gs).grid.cells.centroids(:, 2), N, [])/(milli*meter);

figure('position', [1290, 755, 1275, 559])

val = state.(gs).massfractions{1};
Z = reshape(val, N, []);

surf(X, Y, Z, 'edgecolor', 'none');
colorbar
title('Mass Fraction H2O');
xlabel('x / mm')
ylabel('y / mm')
view([50, 51]);

figure('position', [1290, 755, 1275, 559])

val = state.(gs).massfractions{2};
Z = reshape(val, N, []);

surf(X, Y, Z, 'edgecolor', 'none');
colorbar
title('Mass Fraction O2');
xlabel('x / mm')
ylabel('y / mm')
view([50, 51]);

figure('position', [1290, 755, 1275, 559])

val = state.(gs).pressure;
Z = reshape(val, N, []);

surf(X, Y, Z/barsa, 'edgecolor', 'none');
colorbar
title('Pressure / bar');
xlabel('x / mm')
ylabel('y / mm')
view([50, 51]);

figure('position', [1290, 755, 1275, 559])

val = state.(gs).density;
Z = reshape(val, N, []);

surf(X, Y, Z, 'edgecolor', 'none');
colorbar
title('Density kg/m^3');
xlabel('x / mm')
ylabel('y / mm')
view([50, 51]);


% Current in anode

i = state.Electrolyser.Anode.i;

ind   = model.Electrolyser.couplingTerms{1}.couplingfaces(:, 2);
yc    = model.Electrolyser.Electrolyte.grid.faces.centroids(ind, 2);
areas = model.Electrolyser.Electrolyte.grid.faces.areas(ind);

u = ampere/((centi*meter)^2);
i = (i./areas)/u;

figure
plot(yc/(milli*meter), i);
title('Current in Anode / A/cm^2')
xlabel('height / mm')


% Current in anode

iHp = state.Electrolyser.Anode.iHp;

ind   = model.Electrolyser.couplingTerms{1}.couplingfaces(:, 2);
yc    = model.Electrolyser.Electrolyte.grid.faces.centroids(ind, 2);
areas = model.Electrolyser.Electrolyte.grid.faces.areas(ind);

u = ampere/((centi*meter)^2);
iHp = (iHp./areas)/u;

figure
plot(yc/(milli*meter), iHp);
title('iHp in Anode / A/cm^2')
xlabel('height / mm')

% Faradic effect in Anode

drivingForces.src = @(time) pmControlFunc(time, I, timeswitch, totaltime, 'order', 'I-first');
state = model.evalVarName(state, 'Electrolyser.Anode.iHp', {{'drivingForces', drivingForces}});

i   = state.Electrolyser.Anode.i;
iHp = state.Electrolyser.Anode.iHp;

ind = model.Electrolyser.couplingTerms{1}.couplingfaces(:, 2);
yc  = model.Electrolyser.Electrolyte.grid.faces.centroids(ind, 2);

figure
plot(yc/(milli*meter), iHp./i);
title('Faradic efficiency')
xlabel('height / mm')


%%

figure
plotToolbar(model.GasSupply.grid, states);

%%
