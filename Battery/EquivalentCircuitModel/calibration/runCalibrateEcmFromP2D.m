%% Calibrate ECM model from P2D model
%

soc_step = 0.5;

%%
% The results of the calibration is given in a struct where the fields are functions in battmo table format


jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));

jsonstruct_p2d = mergeJsonStructs({jsonstruct_material, ...
                                   jsonstruct_geometry});
   
jsonstruct = calibrateEcmFromP2D(jsonstruct_p2d, soc_step);

%% 
% plot the results
soc = linspace(0.1, 1.0, 100);
figure
tiledlayout(2, 3, 'tileindexing', 'columnmajor');

nexttile

% the function |setupFunction| creates a function handle from the battmo function input format.
fn = setupFunction(jsonstruct.OCP)
plot(soc, fn(soc))
title('OCP')
xlabel('SOC / -')
ylabel('OCP / V');

nexttile
fn = setupFunction(jsonstruct.R0)
plot(soc, fn(soc))
title('R0')
xlabel('SOC / -')
ylabel('R0 / Ohm');

nexttile
fn = setupFunction(jsonstruct.R1)
plot(soc, fn(soc))
title('R1')
xlabel('SOC / -')
ylabel('R1 / Ohm');

nexttile
fn = setupFunction(jsonstruct.C1)
plot(soc, fn(soc))
title('C1')
xlabel('SOC / -')
ylabel('C1 / Ohm');

nexttile
fn = setupFunction(jsonstruct.R2)
plot(soc, fn(soc))
title('R2')
xlabel('SOC / -')
ylabel('R2 / Ohm');

nexttile
fn = setupFunction(jsonstruct.C2)
plot(soc, fn(soc))
title('C2')
xlabel('SOC / -')
ylabel('C2 / Ohm');




