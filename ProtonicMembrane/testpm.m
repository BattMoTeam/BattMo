clear all
% close all

mrstModule add ad-core

an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';
ctrl  = 'Control';

filename = '/home/xavier/Matlab/Projects/battmo/ProtonicMembrane/protonicMembrane.json';
jsonstruct = fileread(filename);
jsonstruct = jsondecode(jsonstruct);

jsonstruct.(elyte).N = 100;

paramobj = ProtonicMembraneCellInputParams(jsonstruct);

paramobj = setupProtonicMembraneCellGrid(paramobj, jsonstruct);

% Setup model
model = ProtonicMembraneCell(paramobj);

model = model.setupComputationalGraph();

% compute and get computationalGraph (just used for postprocessing)
model = model.validateModel();
cgt   = model.computationalGraph;

model.verbose = true;

% Setup initial state
state0 = model.setupInitialState();

% Setup schedule

tswitch = 1;
T       = 2; % This is not a real time scale, as all the model deals with equilibrium

N1  = 10;
dt1 = tswitch/N1;
N2  = 10;
dt2 = (T - tswitch)/N2;

step.val = [dt1*ones(N1, 1); dt2*ones(N2, 1)];
step.control = ones(numel(step.val), 1);

Imax = 1e1;

control.src   = @(time) controlfunc(time, Imax, tswitch, T, 'order', 'I-first');

schedule = struct('control', control, 'step', step); 

nls = NonLinearSolver();
nls.maxIterations = 20;
nls.errorOnFailure = false;

model.nonlinearTolerance = 1e-7;

[~, states, report] = simulateScheduleAD(state0, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

%%


dothisplot = true;
if dothisplot
    if numel(states) == 0
        return
    end
    
    state = states{end};
    state = model.addVariables(state, control);
    figure(1)
    plot(state.(elyte).pi)
    title('pi')
    figure(2)
    plot(state.(elyte).pi - state.(elyte).phi)
    title('E')
    figure(3)
    plot(state.(elyte).sigmaEl)
    title('sigmaEl')
    figure(4)
    plot(state.(elyte).phi)
    title('phi')

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


