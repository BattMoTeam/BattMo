clear all
close all

% setup mrst modules
mrstModule add ad-core multimodel mrst-gui battery
mrstVerbose off

paramobj = BatteryInputParams();
paramobj = setupBatteryInputParams1D(paramobj);


