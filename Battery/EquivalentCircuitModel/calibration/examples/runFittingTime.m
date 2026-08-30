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
    

