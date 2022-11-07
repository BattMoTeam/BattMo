mrstDebug(20);

mrstModule add ad-core mrst-gui mpfa agmg linearsolvers

% load json struct for the material properties
jsonfilename = fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell', ...
                             'lithium_ion_battery_nmc_graphite.json');

jsonstruct1 = parseBattmoJson(jsonfilename);

% load json struct for the geometrical properties

jsonfilename = fullfile('Examples','utils', 'data', 'geometry1d.json');

jsonstruct2 = parseBattmoJson(jsonfilename);

% merge the two json structs 
jsonstruct = mergeJsonStructs({jsonstruct1, jsonstruct2});

runBatteryJson(jsonstruct);



