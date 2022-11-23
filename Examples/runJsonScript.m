mrstDebug(20);

mrstModule add ad-core mrst-gui mpfa agmg linearsolvers

% load json struct for the material properties
jsonfilename = fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell', ...
                        'lithium_ion_battery_nmc_graphite.json');

jsonstruct = parseBattmoJson(jsonfilename);

% load json struct for the geometrical properties

jsonfilename = fullfile('Examples','utils', 'data', 'geometry1d.json');
jsonstruct_geometry = parseBattmoJson(jsonfilename);

% merge the two json structs 
jsonstruct = mergeJsonStructs({jsonstruct, jsonstruct_geometry});

CRate = jsonstruct.Control.CRate;

jsonstruct.TimeStepping.totalTime = 1.4*hour/CRate;
jsonstruct.TimeStepping.N = 100;

% Run battery simulation with function that takes json input
output = runBatteryJson(jsonstruct);

%%

E             = output.E;
energyDensity = output.energyDensity;
energy        = output.energy;

figure
plot(energyDensity, E)
xlabel('Energy Density [Wh/L]');
ylabel('Voltage [V]');

figure
plot(energy, E)
xlabel('Energy [Wh]');
ylabel('Voltage [V]');


