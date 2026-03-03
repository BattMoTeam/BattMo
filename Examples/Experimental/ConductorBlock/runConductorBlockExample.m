mrstModule add ad-core mrst-gui mpfa

el      = 'ElectronicModel';
thermal = 'ThermalModel';
ctrl    = 'Control';

filename = fullfile(battmoDir(), 'Examples', 'Experimental', 'ConductorBlock', 'conductorBlock.json');
jsonstruct = parseBattmoJson(filename);

jsonstruct.(ctrl).Imax = 0.8;

inputparams = ConductorBlockInputParams(jsonstruct);

gen = ConductorBlockGridGenerator();

gen.xlength = 1.1;
gen.ylength = 1.3;
gen.zlength = 1.5;

gen.nx = 5;
gen.ny = 4;
gen.nz = 3;

inputparams = gen.updateGridInputParams(inputparams, jsonstruct);

model = ConductorBlock(inputparams);

initstate = model.setupInitialState();

totalTime = 1*minute;
N = 10;

step = struct('val'    , totalTime/N*ones(N, 1), ...
              'control', ones(N, 1));

control = struct('src', []);

schedule = struct('step'   , step, ...
                  'control', control);

%% Running the simulation

model.verbose = true;
[~, states, report] = simulateScheduleAD(initstate, model, schedule);

%% plotting

close all

state = states{end};
G = model.(el).grid;

figure
plotToolbar(model.(el).grid, state.(el).phi);

figure
plotToolbar(model.(thermal).grid, state.(thermal).T);

state = model.addVariables(state);
figure
plotToolbar(model.(thermal).grid, state.(thermal).jHeatSource);
