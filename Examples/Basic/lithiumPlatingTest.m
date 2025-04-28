clear all

jsonstruct_material = parseBattmoJson('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json');

jsonstruct_geometry = parseBattmoJson('Examples/jsondatafiles/geometry1d.json');

jsonstruct_timestepping = parseBattmoJson("Examples/Documentation/jsonfiles/Example/timeStepping.json");

jsonstruct_control = parseBattmoJson('Examples/Documentation/jsonfiles/Example/control.json');

jsonstruct = mergeJsonStructs({jsonstruct_material    , ...
                               jsonstruct_geometry    , ...
                               jsonstruct_timestepping, ...
                               jsonstruct_control});


jsonstruct.TimeStepping.totalTime         = 3960/jsonstruct.Control.DRate;
jsonstruct.TimeStepping.numberOfTimeSteps = jsonstruct.TimeStepping.totalTime/100;
jsonstruct.Control.DRate                  = 0.1;

ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
co      = 'Coating';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
am      = 'ActiveMaterial';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
cc      = 'CurrentCollector';


jsonstruct.(ne).(co).(am).useLithiumPlating = true;

[model, inputparams, jsonstruct, gridGenerator] = setupModelFromJson(jsonstruct);





