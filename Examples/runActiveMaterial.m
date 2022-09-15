%% run stand-alone active material model

% clear the workspace and close open figures
clear
close all

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa

%% Setup the properties of Li-ion battery materials and cell design
jsonstruct = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));

paramobj = BatteryInputParams(jsonstruct);

% Some shorthands used for the sub-models
ne    = 'NegativeElectrode';
pe    = 'PositiveElectrode';
am    = 'ActiveMaterial';
sd    = 'SolidDiffusion';
itf   = 'Interface';
elyte = 'Electrolyte';

%% We setup the battery geometry ("bare" battery with no current collector). We use that to recover the parameters for the active material of the positive electrode, which is instantiate later
gen = BareBatteryGenerator3D();
% We update pamobj with grid data
paramobj = gen.updateBatteryInputParams(paramobj);

paramobj.(pe).(am).InterDiffusionCoefficient = 0;
paramobj.(pe).(am).diffusionModelType = 'full';

paramobj = paramobj.(pe).(am);


paramobj.externalCouplingTerm = [];
paramobj.(sd).np = 1;
xlength = 57e-6; 
G = cartGrid(1, xlength);
G = computeGeometry(G);
paramobj.G = G;

model = ActiveMaterial(paramobj);

dograph = true;

if dograph
    model.isRoot = true;
    cgf = ComputationalGraphFilter(model);
    inspectGraphScript(model);
end

return
%% Setup initial state

N = model.(sd).N;

cElectrodeInit   = 40*mol/litre;
phiElectrodeInit = 3.5;
cElectrolyte     = 5e-1*mol/litre;
phiElectrolyte   = 0;
T                = 298;

% set primary variables
initState.phi           = phiElectrodeInit;
initState.(sd).c        = cElectrodeInit*ones(N, 1);
initState.(sd).cSurface = cElectrodeInit;

% set static variable fields
initState.T = T;
initState.(itf).cElectrolyte   = cElectrolyte;
initState.(itf).phiElectrolyte = phiElectrolyte;

%% setup schedule

total = 1*hour;
n     = 10;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.src = 1e1;

schedule = struct('control', control, 'step', step); 

%% Run simulation
model.verbose = true;
[wellSols, states, report] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps', true); 


%% plotting

time = cellfun(@(state) state.time, states);
cSurface = cellfun(@(state) state.(sd).cSurface, states);

plot(time, cSurface);
