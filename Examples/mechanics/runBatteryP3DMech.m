%% Example where we combine a battery discharge simulation with mechanical simulation
%  We use a a simple relationship between Lithium concentration and load stress in the electrode

clear
close all


%% Import the required modules from MRST

mrstModule add ad-core mrst-gui mpfa

%% Setup the properties of Li-ion battery materials and cell design

jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
jsonstruct.include_current_collectors = true;

% We define some shorthand names for simplicity.
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
co      = 'Coating';
am      = 'ActiveMaterial';
sd      = 'SolidDiffusion';
cc      = 'CurrentCollector';
elyte   = 'Electrolyte';
sep     = 'Separator';
thermal = 'ThermalModel';
ctrl    = 'Control';

% jsonstruct.(ne).(co).(am).diffusionModelType = 'simple';
% jsonstruct.(pe).(co).(am).diffusionModelType = 'simple';

paramobj = BatteryInputParams(jsonstruct);

%% Setup the geometry and computational grid

gen = BatteryGeneratorP3D();

% We change the properties of the grid and update paramobj with the results grid
gen.ylength = 1e-4;
gen.ny      = 10;
paramobj    = gen.updateBatteryInputParams(paramobj);

paramobj.(ne).(cc).effectiveElectronicConductivity = 1e5;
paramobj.(pe).(cc).effectiveElectronicConductivity = 1e5;

%%  Initialize the battery model.

model = Battery(paramobj);

%% Plot the grid

plotBatteryGrid(model);

%% Compute the nominal cell capacity and choose a C-Rate

C      = computeCellCapacity(model);
CRate  = model.(ctrl).CRate;
inputI = (C/hour)*CRate; % current

%% Setup the time step schedule
% Smaller time steps are used to ramp up the current from zero to its operational value. Larger time steps are then used
% for the normal operation.

n           = 25;
dt          = [];
dt          = [dt; repmat(0.5e-4, n, 1).*1.5.^[1:n]'];
totalTime   = 1.4*hour/CRate;
n           = 40;
dt          = [dt; repmat(totalTime/n, n, 1)];
times       = [0; cumsum(dt)];
tt          = times(2 : end);
step        = struct('val', diff(times), 'control', ones(numel(tt), 1));

%% Setup the operating limits for the cell
% The maximum and minimum voltage limits for the cell are defined using stopping and source functions. A stopping
% function is used to set the lower voltage cutoff limit. A source function is used to set the upper voltage cutoff
% limit.

% we setup the control by assigning a source and stop function.
control = model.(ctrl).setupScheduleControl();

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step);

%% Setup the initial state of the model

initstate = model.setupInitialState();

%% Setup the properties of the nonlinear solver
nls = NonLinearSolver();
% Change default maximum iteration number in nonlinear solver
nls.maxIterations = 10;
% Change default behavior of nonlinear solver, in case of error
nls.errorOnFailure = false;
% Timestep selector
nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', {{ctrl, 'E'}}, ...
                                                  'targetChangeAbs', 0.03);
% Change default tolerance for nonlinear solver
model.nonlinearTolerance = 1e-5;
% Set verbosity of the solver (if true, value of the residuals for every equation is given)
model.verbose = true;

%% Run simulation

[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, ...
                                                'OutputMinisteps', true, ...
                                                'NonLinearSolver', nls);

%%  Process output and recover the output voltage and current from the output states.

ind = cellfun(@(x) not(isempty(x)), states);
states = states(ind);

E    = cellfun(@(x) x.(ctrl).E, states);
I    = cellfun(@(x) x.(ctrl).I, states);
time = cellfun(@(x) x.time, states);

%% Setup of the mechanical problem

mrstModule add vemmech

opt = struct('E',1,'nu',0.3);
G = createAugmentedGrid(model.G);
G = computeGeometry(G);

%% Mechanical properties

Ev  = repmat(opt.E, G.cells.num, 1);
nuv = repmat(opt.nu, G.cells.num, 1);
C   = Enu2C(Ev, nuv, G);

%% We setup the boundary conditions for the mechanical problem

bc  = cell(2*G.griddim, 1);
if(G.griddim == 2)
    oside = {'Left', 'Right', 'Back', 'Front'};
else
    oside = {'Left', 'Right', 'Back', 'Front', 'Bottom', 'Top'};
end

for i = 1:numel(oside);
    bc{i} = pside([], G, oside{i}, 0);
    bc{i} = rmfield(bc{i}, 'type');
    bc{i} = rmfield(bc{i}, 'sat');
end

for i = 1:numel(bc);
    inodes = mcolon(G.faces.nodePos(bc{i}.face), G.faces.nodePos(bc{i}.face+1)-1);
    nodes = unique(G.faces.nodes(inodes));
    bc{i}.el_bc = struct('disp_bc', struct('nodes', nodes, 'uu', 0, 'faces', bc{i}.face, 'uu_face', 0, 'mask', true(numel(nodes), G.griddim)), ...
        'force_bc', []);
end

bc = {bc{1}, bc{3}, bc{4}};
bc{2}.el_bc.disp_bc.mask(:, 1) = 0;
bc{3}.el_bc.disp_bc.mask(:, 1) = 0;
nodes = [];
faces = [];
mask = [];

for i = 1:numel(bc)
    nodes = [nodes; bc{i}.el_bc.disp_bc.nodes];
    faces = [faces; bc{i}.el_bc.disp_bc.faces];
    mask  = [mask; bc{i}.el_bc.disp_bc.mask];
end

bcdisp =@(x) x*0;
disp_node = bcdisp(G.nodes.coords(nodes, :));
disp_faces = bcdisp(G.faces.centroids(faces, :));
el_bc = struct('disp_bc', struct('nodes'  , nodes     , ...
                                 'uu'     , disp_node , ...
                                 'faces'  , faces     , ...
                                 'uu_face', disp_faces, ...
                                 'mask'   , mask), ...
               'force_bc', []);
load = @(x) -repmat([0, 0], size(x, 1), 1);

%% We run the mechanical simulation for each time step

state0 = initstate;
figure
plotGrid(G)
ax      = axis();
ax(2)   = ax(2)*1.5;
ax(end) = ax(end)*1.1;

state0 = model.evalVarName(state0, {ne, co, am, sd, 'cAverage'});
state0 = model.evalVarName(state0, {pe, co, am, sd, 'cAverage'});

for i = 1 : numel(states)

    state = states{i};

    state = model.evalVarName(state, {ne, co, am, sd, 'cAverage'});
    state = model.evalVarName(state, {pe, co, am, sd, 'cAverage'});

    T = state.ThermalModel.T - state0.ThermalModel.T;

    cna = zeros(G.cells.num, 1);
    cpa = zeros(G.cells.num, 1);

    ind      = model.(ne).(co).G.mappings.cellmap;
    cna(ind) = state.(ne).(co).(am).(sd).cAverage - state0.(ne).(co).(am).(sd).cAverage;
    ind      = model.(pe).(co).G.mappings.cellmap;
    cpa(ind) = state.(pe).(co).(am).(sd).cAverage - state0.(pe).(co).(am).(sd).cAverage;

    dc = (cna + cpa*0.5)/1000;

    [uVEM, extra] = VEM_linElast(G, C, el_bc, load, ...
                                 'pressure', dc*0.05, ...
                                 'experimental_scaling',false);

    vdiv = VEM_div(G);
    mdiv = vdiv*reshape(uVEM', [], 1)./G.cells.volumes;

    clf
    plotCellDataDeformed(G, mdiv, uVEM, 'EdgeAlpha',0.04);
    colorbar();
    for k = 1:numel(el_bc)
        plotFaces2D(G, el_bc.disp_bc.faces);
    end
    axis(ax)
    pause(0.1)

end

%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
