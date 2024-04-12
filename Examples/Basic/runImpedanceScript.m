clear all
close all

mrstModule add ad-core mrst-gui mpfa agmg linearsolvers

jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));

% We define some shorthand names for simplicity.
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
co      = 'Coating';
am      = 'ActiveMaterial';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
cc      = 'CurrentCollector';

jsonstruct.use_thermal = false;
jsonstruct.include_current_collectors = false;

inputparams = BatteryInputParams(jsonstruct);
gen = BatteryGeneratorP2D();

inputparams = gen.updateBatteryInputParams(inputparams);

figure
hold on

socs = linspace(0.1, 1, 5);

for isoc = 1 : numel(socs)

    soc = socs(isoc);
    
    impsolv = ImpedanceSolver(inputparams, soc);

    omegas = linspace(-3, 6, 500);
    omegas = 10.^omegas;
    Z = impsolv.computeImpedance(omegas);
    
    Z = (Z*gen.faceArea)/((centi*meter)^2);

    legtxt = sprintf('SOC=%g', soc);
    plot(real(Z), -imag(Z), 'displayname', legtxt);
    xlabel('real(Z) / Ωcm^2')
    ylabel('-imag(Z) / Ωcm^2')

end

legend show
title('Impedance')
