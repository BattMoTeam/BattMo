%% Setup material properties

jsonfilename = fullfile('Examples', 'Experimental', 'RegionBruggeman', 'region_bruggeman.json');
jsonstruct_bruggeman = parseBattmoJson(jsonfilename);

jsonfilename = fullfile('ParameterData'        , ...
                        'BatteryCellParameters', ...
                        'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_nmc_graphite.json');
jsonstruct_material = parseBattmoJson(jsonfilename);

jsonstruct_material.include_current_collectors = true;

%% Setup geometry
jsonfilename = fullfile('Examples'     , ...
                        'JsonDataFiles', ...
                        'geometry3d.json');
jsonstruct_geometry = parseBattmoJson(jsonfilename);

%% Setup Control
jsonfilename = fullfile('Examples', 'JsonDataFiles', 'cc_discharge_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);

jsonstruct = mergeJsonStructs({jsonstruct_bruggeman, ...
                               jsonstruct_geometry , ...
                               jsonstruct_material , ...
                               jsonstruct_control});


output = runBatteryJson(jsonstruct);
