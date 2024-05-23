
%% prepare input for co2membrane

filename = 'CO2membrane/co2membrane.json';
jsonstruct = parseBattmoJson(filename);

inputparams = CO2membraneInputParams(jsonstruct);

gen = CO2membraneGridGenerator();

inputparams = gen.updateInputParams(inputparams);

%% Setup model for membrane side (using Feed input data)

model = CO2membraneSide(inputparams.Feed);
model.isRootSimulationModel = true;

shortNames =  {'Boundary' ,'bd'  ;
               'bcMolFractionDefinitions', 'bcmfdef';
               'bcFluxDefinition', 'bcfluxdef'};

model = model.equipModelForComputation('shortNames', shortNames);

%%

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
