clear all
close all

mrstModule add ad-core mrst-gui mpfa agmg linearsolvers
addpath('C:\Users\Alexandre Fichter\Documents\stage_3A\contenu stage\data_August\jp3-eis');

filename = 'C:\Users\Alexandre Fichter\Documents\stage_3A\contenu stage\2026-alexandre-fichter\jp3-eis\jp3-opt-1d-full.json';
jsonstruct = parseBattmoJson(filename);

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

set(0, 'defaultlinelinewidth', 3);


params0 = [3e-3, 0.00260, 100, 0.00132, 1];

filename = 'C:\Users\Alexandre Fichter\Documents\stage_3A\contenu stage\data_August\ank_data\Supplementary material\02_Electrical_characterization\EIS\131-828_EIS_01_MB_CD8.txt';
% filename = '/home/xavier/Matlab/Projects/battmo/Data/131-828_EIS_01_MB_CD8.txt';

[Z_re_exp, Z_im_exp, omega] = load_experimental_data(filename);

scales = [10e-2, 10e-3, 5000, 100e-3, 50];

feis = FittingEIS(params0, scales, Z_re_exp, Z_im_exp, omega);

[~, ~, best_params, fitting_error] = feis.optimizationBFGS();

[Z_re_fit, Z_im_fit] = load_nyquist(best_params, feis.omega);

%% Comparison with Ank's experimental data

figure


subplot(3,1,1);
legtxt = sprintf('exp data');
semilogx(feis.omega, feis.Z_re_exp,'r', 'displayname', legtxt);
hold on;
legtxt = sprintf('fitted data');
semilogx(feis.omega, Z_re_fit, 'b', 'displayname', legtxt); 
xlabel('Omega');
ylabel('Z_{re} '); 

subplot(3,1,2);
legtxt = sprintf('exp data');
semilogx(feis.omega, feis.Z_im_exp, 'r',  'displayname', legtxt);
hold on;
legtxt = sprintf('fitted data');
semilogx(feis.omega, Z_im_fit, 'b', 'displayname', legtxt); 
xlabel('Omega');
ylabel('-Z_{im} '); 

subplot(3,1,3);
legtxt = sprintf('exp data');
plot(feis.Z_re_exp, feis.Z_im_exp, 'r',  'displayname', legtxt);
hold on;
legtxt = sprintf('fitted data');
plot(feis.Z_re_exp, Z_im_fit, 'b',  'displayname', legtxt);
axis equal;




hold on

socs = linspace(0.1, 1, 2);

socs = [socs, 0.2];

for isoc = 1 : numel(socs)
    soc = socs(isoc);
    
    impsolv = ImpedanceSolver(inputparams, 'soc', soc, 'computeSteadyState', false);

    omegas = linspace(-2, 4, 500);
    omegas = 10.^omegas;
    Z = impsolv.computeImpedance(omegas);
    
    % Z = (Z*gen.faceArea)/((centi*meter)^2);
    Z = (Z*gen.faceArea);


    subplot(3,1,1);
    legtxt = sprintf('SOC=%g', soc);
    semilogx(omegas, real(Z), 'displayname', legtxt);
    

    subplot(3,1,2);
    legtxt = sprintf('SOC=%g', soc);
    semilogx(omegas, -imag(Z), 'displayname', legtxt);
    
    subplot(3,1,3);
    legtxt = sprintf('SOC=%g', soc);
    plot(real(Z), -imag(Z), 'displayname', legtxt);
    xlabel('real(Z) / Ωcm^2')
    ylabel('-imag(Z) / Ωcm^2')

end

legend show
title('Impedance')
axis equal
