%% Alkaline Membrane Electrolyser

%% Setup input
% Setup the physical properties for the electrolyser using json input file :battmofile:`alkalineElectrolyser.json<Electrolyser/Parameters/alkalineElectrolyser.json>`

jsonstruct_material = parseBattmoJson('Electrolyser/Parameters/alkalineElectrolyser.json');
jsonstruct_geometry = parseBattmoJson('Electrolyser/Parameters/electrolysergeometry1d.json');

jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                               jsonstruct_geometry});

inputparams = ElectrolyserInputParams(jsonstruct);
inputparams = setupElectrolyserGridFromJson(inputparams, jsonstruct);


options = [];
options.stateInitialization.initializationSetup = 'given state';
options.stateInitialization.computeSteadyState  = false;

extrastructs = [];
extrastructs.initstate = initstate;

impsolv = ElectrolyserImpedanceSolver(inputparams, options, extrastructs);






