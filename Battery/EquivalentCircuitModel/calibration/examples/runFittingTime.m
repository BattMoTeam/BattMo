%% Time fitting


%% Setup current input
time_init = (0:10:1000)';
current_init = ones(size(time_init));

current_init(time_init >= 0   & time_init < 200)   = 2.5;  % 1. Décharge (0 à 200s)
current_init(time_init >= 200 & time_init < 400)   = -1.5; % 2. Charge (200 à 400s)
current_init(time_init >= 400 & time_init < 700)   = 0.0;  % 3. Repos (400 à 700s)
current_init(time_init >= 700 & time_init <= 1000) = 3.5;  % 4. Décharge (700 à 1000s)

jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));
jsonstruct = mergeStructs({jsonstruct_material, jsonstruct_geometry});

params0 = [0.05052, 1.12673, 59119.9, 0.03155, 11054.0];  % initial condition: C1>2*C2

a = 1000;

pmin = params0 / a;
pmax = params0 * a;
scales = [pmin, pmax];

ftime = FittingTime(jsonstruct, time_init, current_init, params0, scales);

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

