%% run stand-alone active material model

% clear the workspace and close open figures
clear
close all

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa

%% Setup the properties of Li-ion battery materials and cell design
jsonstruct = parseBattmoJson('ParameterData/ParameterSets/Chen2020/chen2020_lithium_ion_battery.json');

% Some shorthands used for the sub-models
ne  = 'NegativeElectrode';
pe  = 'PositiveElectrode';
am  = 'ActiveMaterial';
sd  = 'SolidDiffusion';
itf = 'Interface';
sei = 'SolidElectrodeInterface';
sr  = 'SideReaction';
elyte = 'Electrolyte';

%% We setup the battery geometry ("bare" battery with no current collector). We use that to recover the parameters for the active material of the positive electrode, which is instantiate later
paramobj = SEIActiveMaterialInputParams(jsonstruct.(pe).(am));

paramobj.(sd).useSimplifiedDiffusionModel = false;
paramobj.InterDiffusionCoefficient = 0;
paramobj.externalCouplingTerm = [];
paramobj.(sd).N   = 10;
paramobj.(sd).np  = 1;
paramobj.(sei).N  = 10;
paramobj.(sei).np = 1;
xlength = 57e-6; 
G = cartGrid(1, xlength);
G = computeGeometry(G);
paramobj.G = G;

model = SEIActiveMaterial(paramobj);


%% Setup initial state

Nsd = model.(sd).N;
Nsei = model.(sei).N;

cElectrodeInit   = 40*mol/litre;
phiElectrodeInit = 3.5;
cElectrolyte     = 5e-1*mol/litre;
phiElectrolyte   = 0;
T                = 298;
cExternal        = 4.541*mol/litre;

% set primary variables
initState.phi              = phiElectrodeInit;
initState.(sd).c           = cElectrodeInit*ones(Nsd, 1);
initState.(sd).cSurface    = cElectrodeInit;
initState.(sei).c          = cElectrolyte*ones(Nsei, 1);
initState.(sei).cInterface = cElectrolyte;
initState.(sei).delta      = 5e-9;
initState.R                = 0;

% set static variable fields
initState.T = T;
initState.(itf).cElectrolyte   = cElectrolyte;
initState.(itf).phiElectrolyte = phiElectrolyte;
initState.(sr).phiElectrolyte  = phiElectrolyte;
initState.(sei).cExternal      = cExternal;


%% setup schedule

total = 1*hour;
n     = 10;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.src = 1e-6;

schedule = struct('control', control, 'step', step); 

%% Run simulation

model.verbose = true;
[wellSols, states, report] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps', true); 



