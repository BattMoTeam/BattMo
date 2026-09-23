%% ECM fitting from time dependent measurements
%


%%
% We use the class FittigTime to fit the parameters of a Thevenin model to given time
% dependent data. The data consists of current and voltage measurement as function of time. Here, we
% generate synthetically this data using a P2D model
%

%% Setup current input
%

clear jsonstruct_control
jsonstruct_control.controlPolicy = "timeControl";

vals = [    0     2.5;
          799   2.5;
          800   -1.5;
         1599   -1.5;
         1600    0;
         2799    0;
         2800    3.5;
         4000   3.5];

clear jsonstruct;
jsonstruct.argumentList = {"time"};
jsonstruct.functionFormat = "tabulated";
jsonstruct.dataX = vals(:, 1);
jsonstruct.dataY = vals(:, 2);

jsonstruct_control.value = jsonstruct;

clear jsonstruct;
jsonstruct.argumentList = {"time"};
jsonstruct.functionFormat = "constant";
jsonstruct.value = 1;

jsonstruct_control.type = jsonstruct;

jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));
jsonstruct = mergeStructs({jsonstruct_material, jsonstruct_geometry});

jsonstruct.Control = jsonstruct_control;
jsonstruct.TimeStepping = struct('totalTime', vals(end, 1), ...
                                 "numberOfTimeSteps", 200);

output = runBattery(jsonstruct);

figure
tiledlayout(1, 2);

nexttile
plot(output.time, output.E)
title('Voltage')
xlabel('time / s')
ylabel('E / V');

nexttile
plot(output.time, output.I)
title('Current')
xlabel('time / s')
ylabel('I / V');

%%
%

params0 = [5e-2; 1e-1; 6e5; 3e-3; 1e5];  

clear inputs
inputs.jsonstruct = jsonstruct;
inputs.output     = output;

clear options
options.useP2Dmodel = true;
options.useSimulationOutput = true;
options.scaleFactor = 1e2;
ftime = FittingTime(params0, inputs, options);


[~, ~, best_params, fitting_error] = ftime.optimizationBFGS();

p = best_params;

%% Printing results

ftime.plotresults_thevenin(best_params, fitting_error);
ftime.printResults(best_params, fitting_error);          
    

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
