clear all
close all

% setup mrst modules
mrstModule add ad-core multimodel mrst-gui battery
mrstVerbose off

paramobj = LithiumBatteryInputParams();

paramobj = setupBatteryInputParams1D(paramobj);

battery = Battery2(paramobj);

