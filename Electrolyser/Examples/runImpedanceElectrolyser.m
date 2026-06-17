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

options.stateInitialization.initializationSetup = 'given current';
options.stateInitialization.computeSteadyState  = true;
options.stateInitialization.I                   = -3*ampere/(centi*meter)^2;

options.stateInitialization.dt = 1;

impsolv = ElectrolyserImpedanceSolver(inputparams, options);

set(0, 'defaultlinelinewidth', 3);

omegas = linspace(-4, 2, 30);
omegas = 10.^omegas;
Z = impsolv.computeImpedance(omegas);

figure
plot(real(Z), -imag(Z));
xlabel('real(Z) / Ω')
ylabel('-imag(Z) / Ω')
legend show
title('Impedance')





