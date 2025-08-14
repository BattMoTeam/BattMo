%% Test the json inputs registered in this directory

filename = 'Examples/Documentation/jsonfiles/Example/battery.json';
jsonstruct_material = parseBattmoJson(filename);

filename = 'Examples/Documentation/jsonfiles/Example/geometry.json';
jsonstruct_geometry = parseBattmoJson(filename);

filename = 'Examples/Documentation/jsonfiles/Example/control.json';
jsonstruct_control = parseBattmoJson(filename);

jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                               jsonstruct_geometry, ...
                               jsonstruct_control});

jsonstruct.include_current_collectors = true;

output = runBatteryJson(jsonstruct);
