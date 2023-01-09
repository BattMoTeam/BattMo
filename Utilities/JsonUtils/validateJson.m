clear all

% Reload the module if changed
reload = false;

if reload
    clear classes
    mod = py.importlib.import_module('validationJsonScript');
    py.importlib.reload(mod);
end

% Json files to validate
jsonfiles = {'ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json',
             'ParameterData/ParameterSets/Xu2015/lfp.json',
             'ParameterData/ParameterSets/Chen2020/chen2020_lithium_ion_battery.json'};

for k = 1:numel(jsonfiles)
    jsonfile = jsonfiles{k};
    is_valid = py.validationJsonScript.validate(jsonfile);
    assert(is_valid)
end
