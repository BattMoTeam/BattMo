%% Alkaline Membrane Electrolyser Impedance computation
%

%% Setup model
% 
jsonstruct_material = parseBattmoJson('Electrolyser/Parameters/alkalineElectrolyser.json');
jsonstruct_geometry = parseBattmoJson('Electrolyser/Parameters/electrolysergeometry1d.json');

jsonstruct = mergeStructs({jsonstruct_material, ...
                               jsonstruct_geometry});

%%
% We setup the parameter structure inputparams for an electolyser
inputparams = ElectrolyserInputParams(jsonstruct);
% We setup the discretization grid. The parameters have been loaded in the json structure
inputparams = setupElectrolyserGridFromJson(inputparams, jsonstruct);

%% Setup Impedance solver 
%

%%
% We use a given current. The steady state is computed at setup.
%

options = [];
options.stateInitialization.initializationSetup = 'given current';
% We use a current density of -3 A/cm^2
options.stateInitialization.I = -3*ampere/(centi*meter)^2;
% The steady state is computed at setup.
options.stateInitialization.computeSteadyState  = true;

impsolv = ElectrolyserImpedanceSolver(inputparams, options);

%% Compute the impedance for a range of frequencies
%
omegas = linspace(-4, 2, 100);
omegas = 10.^omegas;
Z = impsolv.computeImpedance(omegas);

%% Plots of the results
%

set(0, 'defaultlinelinewidth', 3);
figure
plot(real(Z), -imag(Z));
xlabel('real(Z) / Ω')
ylabel('-imag(Z) / Ω')
legend show
title('Impedance')







