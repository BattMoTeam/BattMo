%% Pseudo-Two-Dimensional (P2D) Lithium-Ion Battery Model
% This example demonstrates how to setup a P2D model of a Li-ion battery
% and run a simple simulation.

% Clear the workspace and close open figures
close all

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

jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));

% We define some shorthand names for simplicity.
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

jsonstruct.use_thermal = false;
jsonstruct.include_current_collectors = false;

control = struct('controlPolicy', 'Impedance');
jsonstruct.Control = control;

inputparams = BatteryInputParams(jsonstruct);
gen = BatteryGeneratorP2D();

inputparams = gen.updateBatteryInputParams(inputparams);

model = ImpedanceBattery(inputparams);
cgt = model.cgt;

dopreparestate = true;

if dopreparestate
    state = prepareState(inputparams);
    return
end

omegas = linspace(-2, 4, 50);
omegas = 10.^omegas;

clear Z
for iomega = 1 : numel(omegas)
    omega = omegas(iomega);
    Z(iomega) = model.computeImpedance(state, omega);
end

figure
Z = (Z*gen.faceArea)/((centi*meter)^2);
plot(real(Z), -imag(Z));


function state = prepareState(inputparams, gen)
    
    control = struct('controlPolicy'     , 'CCDischarge', ...
                     'rampupTime'        , 0.1          , ...
                     'DRate'             , 1            , ...
                     'lowerCutoffVoltage', 2.4          , ...
                     'upperCutoffVoltage', 4.1);
    inputparams.Control = CCDischargeControlModelInputParams(control);
    model = Battery(inputparams);
    model.Control.Imax = 0;

    N = 10;
    totalTime = 10*hour;
    dt = rampupTimesteps(totalTime, totalTime/N, 3);

    step.val = dt;
    step.control = ones(numel(dt), 1);

    control = model.Control.setupScheduleControl();

    schedule = struct('step'   , step, ...
                       'control', control);
    
    initstate = model.setupInitialState();
    
    [~, states, report] = simulateScheduleAD(initstate, model, schedule);

    state = states{end};
    
end
