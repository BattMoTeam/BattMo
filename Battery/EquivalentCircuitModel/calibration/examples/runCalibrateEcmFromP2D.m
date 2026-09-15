%% Calibration ECM model from P2D model
%
% We calibrate the ECM model parameters from the EIS response of a P2D model. The calibration is
% done for a set of state of charge. 
%

%% P2D reference model 
% We load the data of the P2D model. The input is the same as for the P2D simulation
%
%%
% We load first the material parameters
filename = fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json');
jsonstruct_material = parseBattmoJson(filename);
%%
% Then, the geometrical parameters
filename = fullfile('Examples', 'JsonDataFiles', 'geometryChen.json');
jsonstruct_geometry = parseBattmoJson(filename);
%%
% We merge the two
jsonstruct_p2d = mergeStructs({jsonstruct_material, ...
                               jsonstruct_geometry});

%%
% The calibration is done at different state of charge. To dei

soc_values = linspace(0, 1, 5);

%%
% The results of the calibration is given in a struct where the fields are functions given in the
% battmo table format (see below)
   
jsonstruct = calibrateEcmFromP2D(jsonstruct_p2d, soc_values);

%% 
% plot the results
soc = sort(soc_values);

figure
tiledlayout(2, 3, 'tileindexing', 'columnmajor');

nexttile

% the function |setupFunction| creates a function handle from the battmo function input format.
fn = setupFunction(jsonstruct.OCP);
plot(soc, fn(soc), '*-')
title('OCP')
xlabel('SOC / -')
ylabel('OCP / V');

nexttile
fn = setupFunction(jsonstruct.R0);
plot(soc, fn(soc), '*-')
title('R0')
xlabel('SOC / -')
ylabel('R0 / \Omega');

nexttile
fn = setupFunction(jsonstruct.R1);
plot(soc, fn(soc), '*-')
title('R1')
xlabel('SOC / -')
ylabel('R1 / \Omega');

nexttile
fn = setupFunction(jsonstruct.C1);
plot(soc, fn(soc), '*-')
title('C1')
xlabel('SOC / -')
ylabel('C1 / F');

nexttile
fn = setupFunction(jsonstruct.R2);
plot(soc, fn(soc), '*-')
title('R2')
xlabel('SOC / -')
ylabel('R2 / \Omega');

nexttile
fn = setupFunction(jsonstruct.C2);
plot(soc, fn(soc), '*-')
title('C2')
xlabel('SOC / -')
ylabel('C2 / F');

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
