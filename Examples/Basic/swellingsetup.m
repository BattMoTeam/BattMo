%% Pseudo-Two-Dimensional (P2D) Lithium-Ion Battery Model
% This example demonstrates how to setup a P2D model of a Li-ion battery
% and run a simple simulation.

% Clear the workspace and close open figures
clear all
close all
clc


%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa agmg linearsolvers

%% Setup the properties of Li-ion battery materials and cell design
% The properties and parameters of the battery cell, including the
% architecture and materials, are set using an instance of
% :class:`BatteryInputParams <Battery.BatteryInputParams>`. This class is
% used to initialize the simulation and it propagates all the parameters
% throughout the submodels. The input parameters can be set manually or
% provided in json format. All the parameters for the model are stored in
% the inputparams object.

jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_silicon.json'));

% We define some shorthand names for simplicity.
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
am      = 'ActiveMaterial';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
cc      = 'CurrentCollector';

jsonstruct.use_thermal = false;
jsonstruct.include_current_collectors = false;

jsonstruct.use_particle_diffusion = true;

inputparams = BatteryInputParams(jsonstruct);

inputparams.(ne).(am).InterDiffusionCoefficient = 0;
inputparams.(pe).(am).InterDiffusionCoefficient = 0;

inputparams.(ne).(am).(sd).N = 5;
inputparams.(pe).(am).(sd).N = 5;

inputparams = inputparams.validateInputParams();

gen = BatteryGenerator1D();

% Now, we update the inputparams with the properties of the mesh. 
inputparams = gen.updateBatteryInputParams(inputparams);

inputparams = inputparams.(ne).(am);

model = SwellingMaterial(inputparams);

model = model.setupComputationalGraph();

cgt = model.computationalGraph;

g = cgt.getComputationalGraph

plot(g)
