clear all
close all

mrstModule add ad-core

an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';
ctrl  = 'Control';

filename = '/home/xavier/Matlab/Projects/battmo/OxideElectrolyte/oxidemembrane.json';
jsonstruct = fileread(filename);
jsonstruct = jsondecode(jsonstruct);

jsonstruct.(elyte).N = 100;

inputparams = OxideMembraneCellInputParams(jsonstruct);

inputparams = setupProtonicMembraneCellGrid(inputparams, jsonstruct);

% Setup model
model = OxideMembraneCell(inputparams);

model = model.setupComputationalGraph();

% compute and get computationalGraph (just used for postprocessing)
model = model.validateModel();
cgit   = model.computationalGraph;

model.verbose = true;

% Setup initial state
state0 = model.setupInitialState();

% Setup schedule


T = 1; % This is not a real time scale, as all the model deals with equilibrium
N = 20;

dt = T/N;

step.val = dt*ones(N, 1);
step.control = ones(numel(step.val), 1);

dolineartest = false;
if dolineartest
    step.val = 0;
    step.control = ones(numel(step.val), 1);
end

% Imax = 1e-2*ampere/((centi*meter)^2);
Imax = 0;

control.src = @(time) controlfunc(time, Imax, [], T, 'order', 'alpha-only');

schedule = struct('control', control, 'step', step); 

nls = NonLinearSolver();
nls.maxIterations = 20;
nls.errorOnFailure = false;

model.nonlinearTolerance = 1e-8;

[~, states, report] = simulateScheduleAD(state0, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

%%

set(0, 'defaultlinelinewidth', 3);
set(0, 'defaultaxesfontsize', 15);

xc = model.(elyte).G.cells.centroids;

dothisplot = true;
if dothisplot
    if numel(states) == 0
        return
    end
    
    state = states{end};
    state = model.addVariables(state, control);
    
    figure(1)
    plot(xc, state.(elyte).ce)
    title('ce')
    xlabel('x [m]')
    figure(2)
    plot(xc, state.(elyte).phi)
    title('phi')
    xlabel('x [m]')
    figure(3)
    plot(xc, state.(elyte).pi)
    title('pi')
    xlabel('x [m]')
    figure(4)
    plot(xc, state.(elyte).gradConcCoefs{1})
    title('dch')
    xlabel('x [m]')    
    figure(5)
    plot(xc, state.(elyte).gradConcCoefs{2})
    title('dce')
    xlabel('x [m]')
    figure(6)
    plot(xc, state.(elyte).ch)
    title('ch')
    xlabel('x [m]')
    figure(3)
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


