%% Setup material properties

jsonfilename = fullfile('Examples', 'Experimental', 'RegionBruggeman', 'region_bruggeman.json');
jsonstruct_bruggeman = parseBattmoJson(jsonfilename);

jsonfilename = fullfile('ParameterData'        , ...
                        'BatteryCellParameters', ...
                        'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_nmc_graphite.json');
jsonstruct_material = parseBattmoJson(jsonfilename);

jsonstruct_material.include_current_collectors = true;
jsonstruct_material.use_thermal = false;

%% Setup geometry
jsonfilename = fullfile('Examples'     , ...
                        'JsonDataFiles', ...
                        'geometry1d.json');
jsonstruct_geometry = parseBattmoJson(jsonfilename);

%% Setup Control
jsonfilename = fullfile('Examples', 'JsonDataFiles', 'cc_discharge_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);

jsonstruct_timestepping.TimeStepping.numberOfTimeSteps = 40;

jsonstruct = mergeJsonStructs({jsonstruct_bruggeman, ...
                               jsonstruct_geometry , ...
                               jsonstruct_material , ...
                               jsonstruct_control  , ...
                               jsonstruct_timestepping});

output = runBatteryJson(jsonstruct);

doplot = true;

if doplot

    states = output.states;
    E = cellfun(@(state) state.Control.E, states);
    t = cellfun(@(state) state.time, states);

    figure
    plot(t, E)
    
end
