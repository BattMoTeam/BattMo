clear all
% close all

mrstModule add ad-core mrst-gui mpfa agmg linearsolvers

jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));

jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                               jsonstruct_geometry});

[model, inputparams, ~, gen] = setupModelFromJson(jsonstruct);

c_ne = 29.866*mol/litre; % initial concentration at negative electrode
c_pe = 17.038*mol/litre; % initial concentration at positive electrode

initstate = initStateChen2020(model, c_ne, c_pe);


impsolv = ImpedanceSolver(inputparams, 'initstate', initstate, 'computeSteadyState', false);

%%

set(0, 'defaultlinelinewidth', 3);

omegas = linspace(-4, 2, 30);
omegas = 10.^omegas;
Z = impsolv.computeImpedance(omegas);

figure
hold on

plot(real(Z), -imag(Z), 'displayname', 'battmo');

docompare = true;
if docompare
    p = fileparts(mfilename('fullpath'));
    data = load(fullfile(p, 'utils', 'pybamm_chen_impedances.mat'));
    Zpybamm = data.impedances;
    plot(real(Zpybamm), -imag(Zpybamm), 'displayname', 'pybamm');
end

xlabel('real(Z) / Ω')
ylabel('-imag(Z) / Ω')
legend show
title('Impedance')


