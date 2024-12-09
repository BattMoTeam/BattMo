%% Protonic Membrane model


%% Setup input
filename = fullfile(battmoDir(), 'ProtonicMembrane', 'protonicMembrane.json');
jsonstruct_material = parseBattmoJson(filename);

filename = fullfile(battmoDir(), 'ProtonicMembrane', '1d-PM-geometry.json');
jsonstruct_geometry = parseBattmoJson(filename);

jsonstruct = mergeJsonStructs({jsonstruct_material, jsonstruct_geometry});


%% Input structure setup
% We setup the input parameter structur

inputparams = ProtonicMembraneCellInputParams(jsonstruct);

%% 
% The json structure has been also updated with some default values and stored in the inputparams variable. We retrieve
% it
jsonstruct = inputparams.jsonstruct;

% We setup the grid, which is done by calling the function :battmo:`setupProtonicMembraneCellGrid`
[inputparams, gen] = setupProtonicMembraneCellGrid(inputparams, jsonstruct);

%% Model setup
% We instantiate the model for the protonic membrane cell
model = ProtonicMembraneCell(inputparams);

%%
% The model is equipped for simulation using the following command (this step may become unnecessary in future versions)
model = model.setupForSimulation();

%% Initial state setup
% We setup the initial state using a default setup included in the model
state0 = model.setupInitialState();

%% Schedule schedule
% we setup schedule, which means the timesteps and also the control we want 

schedule = model.Control.setupSchedule(jsonstruct);

%% Nonlinear solver options
% We have added some extra options to the nonlinear solver
nls = NonLinearSolver();
nls.maxIterations  = 20;
nls.errorOnFailure = false;
nls.verbose        = true;

%%
% We change the default tolerance
model.nonlinearTolerance = 1e-8;

%% Run
% We run the simulation
[~, states, report] = simulateScheduleAD(state0, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

%% Plotting
%

an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';
ctrl  = 'Control';

set(0, 'defaultlinelinewidth', 3);
set(0, 'defaultaxesfontsize', 15);

do1D = true;
if do1D
    N = gen.N;
else
    N = gen.Nx;
end

xc = model.(elyte).grid.cells.centroids(1 : N, 1);

dothisplot = true;
if dothisplot
    if numel(states) == 0
        return
    end
    
    state = states{end};
    state = model.addVariables(state, control);
    figure(1)
    plot(xc, state.(elyte).pi(1 : N))
    title('pi')
    xlabel('x [m]')
    figure(2)
    hold on
    plot(xc, state.(elyte).pi(1 : N) - state.(elyte).phi(1 : N))
    title('E')
    xlabel('x [m]')
    figure(3)
    plot(xc, log(state.(elyte).sigmaEl(1 : N)))
    title('log(sigmaEl)')
    xlabel('x [m]')
    figure(4)
    plot(xc, state.(elyte).phi(1 : N))
    title('phi')
    xlabel('x [m]')
    
    fprintf('min E : %g\nmin sigmaEl : %g\n', ...
            min(state.(elyte).pi - state.(elyte).phi), ...
            min(state.(elyte).sigmaEl));

    return
end

%%

close all
figure
hold on
for istate = 1 : numel(states)
    state = states{istate};
    state = model.addVariables(state, control);
    plot(state.(elyte).E, '*-')
    title('E')
end

figure
hold on
for istate = 1 : numel(states)
    state = states{istate};
    state = model.addVariables(state, control);
    plot(state.(elyte).sigmaEl, '*-')
    title('sigmaEl')
end


return

%%

for istate =  numel(states)
    if istate == 0
        break
    end
    figure
    state = states{istate};
    state = model.addVariables(state, control.src);
    plot(state.(elyte).sigmaEl, '*-')
    title('sigmaEl')
end


