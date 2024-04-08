clear all
close all

mrstModule add ad-core

an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';
ctrl  = 'Control';

filename = '/home/xavier/Matlab/Projects/battmo/ProtonicMembrane/protonicMembrane.json';
jsonstruct = parseBattmoJson(filename);

jsonstruct.(elyte).Nx = 1000;

do1D = true;
if do1D
    jsonstruct.Geometry.type = '1D';
    jsonstruct.(elyte).faceArea = 1;
end

inputparams = ProtonicMembraneCellInputParams(jsonstruct);

[inputparams, gen] = setupProtonicMembraneCellGrid(inputparams, jsonstruct);

% Setup model
model = ProtonicMembraneCell(inputparams);

model = model.setupForSimulation();

cgt   = model.computationalGraph;

model.verbose = true;

% Setup initial state
state0 = model.setupInitialState();

% Setup schedule

tswitch = 1;
T       = 2; % This is not a real time scale, as all the model deals with equilibrium

N1  = 20;
dt1 = tswitch/N1;
N2  = 20;
dt2 = (T - tswitch)/N2;

step.val = [dt1*ones(N1, 1); dt2*ones(N2, 1)];
step.control = ones(numel(step.val), 1);

% Imax = 1e-2*ampere/((centi*meter)^2);
% Imax = 0.3*ampere/((centi*meter)^2);
% Imax = 0.01*ampere/((centi*meter)^2);
Imax = 0;

control.src = @(time) controlfunc(time, Imax, tswitch, T, 'order', 'I-first');

schedule = struct('control', control, 'step', step); 

nls = NonLinearSolver();
nls.maxIterations = 20;
nls.errorOnFailure = false;
nls.verbose = true;

model.nonlinearTolerance = 1e-8;

[~, states, report] = simulateScheduleAD(state0, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

%%

set(0, 'defaultlinelinewidth', 3);
set(0, 'defaultaxesfontsize', 15);

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

