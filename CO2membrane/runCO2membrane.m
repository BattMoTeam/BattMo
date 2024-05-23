filename = 'CO2membrane/co2membrane.json';
jsonstruct = parseBattmoJson(filename);

inputparams = CO2membraneInputParams(jsonstruct);

gen = CO2membraneGridGenerator();

inputparams = gen.updateInputParams(inputparams);

model = CO2membrane(inputparams);



shortNames =  {'Feed' ,'fd'  ;
               'Permeate' ,'pt'  ;
               'Boundary' ,'bd'  ;
               'massConses' ,'mcs'  ;
               'bcMolFractionDefinitions', 'bcmfdef';
               'bcFluxDefinition', 'bcfluxdef'};


model = model.equipModelForComputation('shortNames', shortNames);

initstate = model.setupInitialState();

step = struct('val', 1, ...
              'control', 1);

control = struct('src', []);

schedule = struct('step', step, ...
                  'control', control);
model.verbose = true;

[~, states] = simulateScheduleAD(initstate, model, schedule);

state = states{end};
state = model.addVariables(state);
return
cgp = model.cgp
cgt = model.cgt

close all
figure
cgp.plot;

cgt.printRootVariables;
cgt.printTailVariables;
