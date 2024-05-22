mrstDebug(20);

filename = 'CO2membrane/co2membrane.json';
jsonstruct = parseBattmoJson(filename);

inputparams = CO2membraneInputParams(jsonstruct);

gen = CO2membraneGridGenerator();

inputparams = gen.updateInputParams(inputparams);

model = CO2membraneSide(inputparams.Feed);
model.isRootSimulationModel = true;
model = model.equipModelForComputation;

cgp = model.cgp;
cgt = model.cgt;

% cgt.printRootVariables;
% cgt.printTailVariables;

initstate = model.setupInitialState();

step = struct('val', 1, ...
              'control', 1);

control = struct('src', []);

schedule = struct('step', step, ...
                  'control', control);

model.verbose = true;

[~, states] = simulateScheduleAD(initstate, model, schedule);

state = states{end};
