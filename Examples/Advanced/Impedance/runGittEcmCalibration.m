%% Calibrate directional 2-RC equivalent-circuit models from GITT data
% This example loads stand-alone CSV measurements, estimates directional
% OCV curves from current-off relaxation, fits local 2-RC pulse models, and
% assembles separate charge and discharge BattMo ECM input structures.

clear all
close all

dataFile = fullfile(exampleDirectory, 'data', 'gitt', ...
                    'gitt_measurements.csv');

runFullCalibration = false;
if ~runFullCalibration
    % A deterministic representative fit suitable for an interactive run.
    options.maximumOcvIterations = 1;
    options.maximumOptimizerIterations = 15;
    options.multiStarts = 1;
    options.maximumFitSamples = 200;
end

gittCalibration = calibrateEcmFromGitt(dataFile, options);

disp('Accepted local 2-RC fits:');
disp(gittCalibration.finalFits(:, ...
    {'Direction', 'SOC', 'R0_Ohm', 'R1_Ohm', 'Tau1_s', ...
     'R2_Ohm', 'Tau2_s', 'RMSE_V'}));

disp('Held-out, non-uniform-time-weighted validation scores:');
disp(gittCalibration.validation);

plotGittEcmCalibration(gittCalibration);

% Ready-to-use BattMo configuration structures are returned separately:
chargeEcmJsonstruct = gittCalibration.charge.jsonstruct;
dischargeEcmJsonstruct = gittCalibration.discharge.jsonstruct;
