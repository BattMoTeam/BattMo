%% Time fitting
%

%% Setup current input
time_init = (0:10:1000)';
current_init = ones(size(time_init));

clear jsonstruct_control
jsonstruct_control.controlPolicy = "timeControl";

vals = [0 2.5;
        199  2.5;
        200  -1.5;
        399  -1.5;
        400   0;
        699   0;
        700   3.5;
        1000  3.5];

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
jsonstruct.TimeStepping = struct('totalTime', 1000, ...
                                 "numberOfTimeSteps", 200);

% output = runBattery(jsonstruct);

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

params0 = [0.05052; 1.12673; 59119.9; 0.03155; 11054.0];  % initial condition: C1>2*C2

clear inputs
inputs.jsonstruct = jsonstruct;
inputs.output     = output;

clear options
options.useP2Dmodel = true;
options.useSimulationOutput = true;

ftime = FittingTime(params0, inputs, options);

return

[~, ~, best_params, fitting_error] = ftime.optimizationBFGS();

p = best_params;


%% Printing results
doprint = true;
if doprint
    
    ftime.plotresults_thevenin(best_params, fitting_error);
    
    ftime.printResults(best_params, fitting_error);          

    % Define the 3-pulse test current profile
    time_vec = ftime.time_vec(:); 
    current_test = ones(size(time_vec));
    current_test(time_vec >= 0   & time_vec < 200)   = -5;  
    current_test(time_vec >= 200 & time_vec < 400)   = 5;
    current_test(time_vec >= 400 & time_vec < 700)   = -5;  
    current_test(time_vec >= 700 & time_vec <= 1000) = 0;
    
    ftime.plottest(best_params, current_test)
    
end

