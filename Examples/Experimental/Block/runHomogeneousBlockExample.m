el      = 'ElectronicModel';
thermal = 'ThermalModel';
ctrl    = 'Control';
            
jsonstruct = parseBattmoJson('/home/xavier/Matlab/Projects/battmo/Examples/Experimental/Block/block.json');

jsonstruct.(ctrl).Imax = 0.8;

inputparams = HomogeneousBlockInputParams(jsonstruct);

gen = HomogeneousBlockGridGenerator();

gen.xlength = 1.1;
gen.ylength = 1.3;
gen.zlength = 1.5;

gen.nx = 5;
gen.ny = 7;
gen.nz = 9;

inputparams = gen.updateGridInputParams(inputparams);

model = HomogeneousBlock(inputparams);

cgt = model.cgt;
cgp = model.cgp;

initstate = model.setupInitialState();

totalTime = 1*minute;
N = 10;

step = struct('val'    , totalTime/N*ones(N, 1), ...
              'control', ones(N, 1));

control = struct('src', []);

schedule = struct('step'   , step, ...
                  'control', control);

cgt = model.cgt;
cgp = model.cgp;

%% Running the simulation

model.verbose = true;
[~, states, report] = simulateScheduleAD(initstate, model, schedule);

%% plotting

state = states{end};
G = model.(el).grid;

figure
plotCellData(G, state.(el).phi);

figure
plotToolbar(model.(thermal).grid, state.(thermal).T);
