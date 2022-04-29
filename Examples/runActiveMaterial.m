%% Pseudo-Two-Dimensional (P2D) Lithium-Ion Battery Model
% This example demonstrates how to setup a P2D model of a Li-ion battery
% and run a simple simulation.

% clear the workspace and close open figures
clear
close all
clc


%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa

%% Setup the properties of Li-ion battery materials and cell design
% The properties and parameters of the battery cell, including the
% architecture and materials, are set using an instance of
% :class:`BatteryInputParams <Battery.BatteryInputParams>`. This class is
% used to initialize the simulation and it propagates all the parameters
% throughout the submodels. The input parameters can be set manually or
% provided in json format. All the parameters for the model are stored in
% the paramobj object.
% The input parameters can be given in json format. The json file is read and used to populate the paramobj object.
jsonstruct = parseBattmoJson('ParameterData/ParameterSets/Chen2020/chen2020_lithium_ion_battery.json');

paramobj = BatteryInputParams(jsonstruct);

% Some shorthands used for the sub-models
ne    = 'NegativeElectrode';
pe    = 'PositiveElectrode';
am    = 'ActiveMaterial';
sd    = 'SolidDiffusion';
itf   = 'Interface';
elyte = 'Electrolyte';

%% We setup the battery geometry ("bare" battery with no current collector).
gen = BareBatteryGenerator3D();
% We update pamobj with grid data
paramobj = gen.updateBatteryInputParams(paramobj);

paramobj.(ne).(am).InterDiffusionCoefficient = 0;
paramobj.(pe).(am).InterDiffusionCoefficient = 0;

paramobj.(ne).(am).(sd).useSimplifiedDiffusionModel = false;
paramobj.(pe).(am).(sd).useSimplifiedDiffusionModel = false;

paramobj = paramobj.(pe).(am);

paramobj.externalCouplingTerm = [];
paramobj.(sd).np = 1;
G = cartGrid(1);
G = computeGeometry(G);
paramobj.G = G;

%%  The Battery model is initialized by sending paramobj to the Battery class constructor 

model = ActiveMaterial(paramobj);

N = model.(sd).N;

cElectrodeInit   = 1*mol/litre;
phiElectrodeInit = 0;
cElectrolyte     = 2*mol/litre;
phiElectrolyte   = 1;
T                = 298;

% set primary variables
initState.phi           = phiElectrodeInit;
initState.(sd).c        = cElectrodeInit*ones(N, 1);
initState.(sd).cSurface = cElectrodeInit;

% set static variable fields
initState.T              = T;
initState.jCoupling      = 0;
initState.jFaceCoupling  = 0;
initState.(itf).cElectrolyte   = cElectrolyte;
initState.(itf).phiElectrolyte = phiElectrolyte;

total = 1*minute;
n     = 10;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.src = [];

schedule = struct('control', control, 'step', step); 

[wellSols, states, report] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps', true); 



