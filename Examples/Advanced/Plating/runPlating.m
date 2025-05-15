%% graph plot
jsonfilename = fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_nmc_graphite.json');
jsonstruct = parseBattmoJson(jsonfilename);

ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
co      = 'Coating';
am      = 'ActiveMaterial';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
cc      = 'CurrentCollector';

jsonstruct_am = jsonstruct.(ne).(co).(am);

jsonstruct_am.SEImodel          = 'none';
jsonstruct_am.SolidDiffusion.np = 1;
jsonstruct_am.useLithiumPlating = true;

inputparams = ActiveMaterialInputParams(jsonstruct_am);
model = ActiveMaterial(inputparams);
cgp = model.cgp;
close all
cgp.plot();
