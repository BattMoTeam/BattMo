
%% prepare input for co2membrane

filename = 'CO2capture/CO2capture.json';
jsonstruct = parseBattmoJson(filename);

inputparams = CO2captureInputParams(jsonstruct);

gen = CO2captureGridGenerator();

inputparams = gen.updateInputParams(inputparams);

%% Setup model for membrane channel (using Feed input data)

model = CO2captureChannel(inputparams.Feed);
model.isRootSimulationModel = true;

shortNames =  {'Boundary' ,'bd'  ;
               'bcMolFractionDefinitions', 'bcmfdef';
               'bcFluxDefinition', 'bcfluxdef'};

model = model.equipModelForComputation('shortNames', shortNames);

%%

cgit = model.cgit;

return
%%

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
