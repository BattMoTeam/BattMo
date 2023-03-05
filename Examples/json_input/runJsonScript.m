mrstModule add ad-core mrst-gui mpfa agmg linearsolvers

% load json struct for the material properties
jsonfilename = fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_nmc_graphite.json');
jsonstruct_material = parseBattmoJson(jsonfilename);

% load json struct for the geometrical properties
jsonfilename = fullfile('Examples','utils', 'data', 'geometry1d.json');
jsonstruct_geometry = parseBattmoJson(jsonfilename);

% load json struct for the geometrical properties
jsonfilename = fullfile('Examples','utils', 'data', 'ie_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);

% load json struct for the geometrical properties
jsonfilename = fullfile('Examples','utils', 'data', 'simulation_parameters.json');
jsonstruct_simparams = parseBattmoJson(jsonfilename);

% load json struct for output extra specifications
jsonfilename = fullfile('Examples','utils', 'data', 'extra_output.json');
jsonstruct_output = parseBattmoJson(jsonfilename);

% We merge the json structures. The function issues a warning if a parameter is set with different values in the given
% structures. The rule is that the first value takes precedence.
jsonstruct = mergeJsonStructs({jsonstruct_geometry , ...
                               jsonstruct_material , ...
                               jsonstruct_control  , ...
                               jsonstruct_simparams, ...
                               jsonstruct_output   , ...                               
                              });


%% We adjust the total time with respect to the given CRate.

CRate = jsonstruct.Control.CRate;
jsonstruct.TimeStepping.totalTime = 1.4*hour/CRate;
jsonstruct.TimeStepping.N = 40;

%% Run battery simulation with function that takes json input

output = runBatteryJson(jsonstruct);

%% Plotting

states = output.states;

E = cellfun(@(x) x.Control.E, states); 
I = cellfun(@(x) x.Control.I, states);
time = cellfun(@(x) x.time, states); 

figure()
subplot(1,2,1)
plot(time/hour, E)
xlabel('time [hours]')
ylabel('Cell Voltage [V]')

subplot(1,2,2)
plot(time/hour, I)
xlabel('time [hours]')
ylabel('Cell Current [A]')

