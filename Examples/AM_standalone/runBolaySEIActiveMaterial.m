%% Particle simulation with SEI layer growth (Bolay et al 2022)

% clear the workspace and close open figures
clear
close all

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa

ne    = 'NegativeElectrode';
pe    = 'PositiveElectrode';
am    = 'ActiveMaterial';
co    = 'Coating';
sd    = 'SolidDiffusion';
itf   = 'Interface';
sei   = 'SolidElectrodeInterface';
sr    = 'SideReaction';
elyte = 'Electrolyte';

%% Setup the properties of the Li-ion battery materials and of the cell design
jsonfilename = fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_nmc_graphite.json');
jsonstruct_material = parseBattmoJson(jsonfilename);
jsonstruct = jsonstruct_material.(ne).(co).(am);

jsonfilename = fullfile('ParameterData', 'ParameterSets', 'Bolay2022', 'bolay_sei_interface.json');
jsonstruct_bolay = parseBattmoJson(jsonfilename);


jsonstruct.(itf) = mergeJsonStructs({jsonstruct.(itf), ...
                                     jsonstruct_bolay});


jsonstruct.sei_type              = 'bolay';
jsonstruct.(sd).N                = 10;
jsonstruct.isRootSimulationModel = true;

jsonstruct.(sd).referenceDiffusionCoefficient = 1e-14;

rp = jsonstruct.(sd).particleRadius ;
jsonstruct.(itf).volumetricSurfaceAreas  = 3./rp;


inputparams = ActiveMaterialInputParams(jsonstruct);

% We initiate the model
model = ActiveMaterial(inputparams);

model = model.setupForSimulation();



%% Setup initial state

Nsd  = model.(sd).N;

% Initial concentration value at the electrode
cElectrodeInit = 0.75*model.(itf).saturationConcentration;
% Initial value of the potential at the electrode
phiElectrodeInit = 0;
% Initial concentration value in the electrolyte
cElectrolyte = 5e-1*mol/litre;
% Temperature
T = 298.15;

% The following datas come from :cite:`Bolay2022` (supplementary material)
% Length of SEI layer
SEIlength = 10*nano*meter;
% SEI voltage drop
SEIvoltageDrop = 0;

% We compute the OCP from the given data and use it to assign electrical potential in electrolyte
initState.T = T;
initState.(sd).cSurface = cElectrodeInit;
initState = model.evalVarName(initState, {itf, 'OCP'});

OCP = initState.(itf).OCP;
phiElectrolyte = phiElectrodeInit - OCP;

% From the values computed above we set the values of the initial state
initState.E                    = phiElectrodeInit;
initState.(sd).c               = cElectrodeInit*ones(Nsd, 1);
initState.(itf).SEIlength      = SEIlength;
initState.(itf).SEIvoltageDrop = SEIvoltageDrop;

% We set also static variable fields
initState.(itf).cElectrolyte   = cElectrolyte;
initState.(itf).phiElectrolyte = phiElectrolyte;

%% Setup schedule

Imax = 3e-12*ampere;

scalings = {};
coef = Imax/PhysicalConstants.F;
scalings{end + 1} = {{sd, 'massCons'}, coef};
scalings{end + 1} = {{sd, 'solidDiffusionEq'}, coef};
scalings{end + 1} = {{'chargeCons'}, Imax};

L0 = 1*nano*meter;
k  = model.(itf).SEIionicConductivity;
rp = model.(sd).particleRadius;
F  = PhysicalConstants.F;
vsa = model.(sd).volumetricSurfaceArea;

coef = Imax*L0*k/(4*pi/3*(rp)^3*F*vsa);

scalings{end + 1} = {{itf, 'SEIvoltageDropEquation'}, coef};

De = model.(itf).SEIelectronicDiffusionCoefficient;
ce = model.(itf).SEIintersticialConcentration;

coef = De*ce/L0;

scalings{end + 1} = {{itf, 'SEImassCons'}, coef};

model.scalings = scalings;

total = 60*minute;
n     = 200;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

% rampup value for the current function, see rampupSwitchControl
tup = 1e-2*minute;
srcfunc = @(time) rampupControl(time, tup, Imax);

cmin = (model.(itf).guestStoichiometry0)*(model.(itf).saturationConcentration);
control.stopFunction = @(model, state, state0_inner) (state.(sd).cSurface <= cmin);
control.src = srcfunc;

schedule = struct('control', control, 'step', step);

%% Setup non-linear solver

nls = NonLinearSolver();
nls.errorOnFailure = false;

model.nonlinearTolerance = 1e-5;

%% Run simulation

model.verbose = true;
[~, states, report] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

%% Plotting

set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);
set(0, 'defaultlinelinewidth', 3);
set(0, 'defaultfigureposition', [10, 10, 800, 400]);

ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

time = cellfun(@(state) state.time, states);

cSurface = cellfun(@(state) state.(sd).cSurface, states);
figure
plot(time/hour, cSurface/(1/litre));
xlabel('time / h');
ylabel('Surface concentration / mol/L');
title('Surface concentration');

E = cellfun(@(state) state.E, states);
figure
plot(time/hour, E);
xlabel('time / h');
ylabel('Potential / V');
title('Potential');


cmin = cellfun(@(state) min(state.(sd).c), states);
cmax = cellfun(@(state) max(state.(sd).c), states);

for istate = 1 : numel(states)
    states{istate} = model.evalVarName(states{istate}, {sd, 'cAverage'});
end

caver = cellfun(@(state) max(state.(sd).cAverage), states);

figure
hold on
plot(time/hour, cmin /(mol/litre), 'displayname', 'cmin');
plot(time/hour, cmax /(mol/litre), 'displayname', 'cmax');
plot(time/hour, caver/(mol/litre), 'displayname', 'total concentration');
title('Concentration in particle / mol/L')
legend show

c = states{end}.(sd).c;
r = linspace(0, model.(sd).particleRadius, model.(sd).N);

figure
plot(r, c/(mol/litre));
xlabel('radius / m')
ylabel('concentration / mol/L')
title('Particle concentration profile (last time step)')

seilength = cellfun(@(state) state.(itf).SEIlength, states);

figure
plot(time/hour, seilength);
xlabel('time / hour')
ylabel('length / m')
title('SEI layer length')

SEIvoltageDrops = cellfun(@(state) state.Interface.SEIvoltageDrop, states);

figure;
plot(time/hour, SEIvoltageDrops);
xlabel('time / h');
ylabel('Drop /V');
title('SEI drop');
grid on;
