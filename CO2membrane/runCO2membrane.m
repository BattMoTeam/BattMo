filename = 'CO2membrane/co2membrane.json';
jsonstruct = parseBattmoJson(filename);

inputparams = CO2membraneInputParams(jsonstruct);

gen = CO2membraneGridGenerator();

inputparams = gen.updateInputParams(inputparams);

model = CO2membrane(inputparams);

model = model.equipModelForComputation();

initstate = model.setupInitialState();

step = struct('val', 1, ...
              'control', 1);

control = struct('src', []);

schedule = struct('step', step, ...
                  'control', control);

output = simulateScheduleAD(initstate, model, schedule);

cgp = model.cgp
cgt = model.cgt

close all
figure
cgp.plot;

cgt.printRootVariables;
cgt.printTailVariables;
