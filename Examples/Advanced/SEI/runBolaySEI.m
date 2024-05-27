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
jsonstruct = parseBattmoJson(jsonfilename);

jsonfilename = fullfile('ParameterData', 'ParameterSets', 'Bolay2022', 'bolay_sei_interface.json');
jsonstruct_bolay = parseBattmoJson(jsonfilename);


jsonstruct.(ne).(co).(am) = mergeJsonStructs({jsonstruct.(ne).(co).(am), ...
                                              jsonstruct_bolay});

jsonstruct.(ne).(co).(am).SEImodel                     = 'Bolay';

inputparams = BatteryInputParams(jsonstruct);

gen = BatteryGeneratorP2D();
inputparams = gen.updateBatteryInputParams(inputparams);


model = GenericBattery(inputparams);
model = model.equipModelForComputation();

%% Setup the schedule
%

timestep.numberOfTimeSteps = 100;

step    = model.Control.setupScheduleStep(timestep);
control = model.Control.setupScheduleControl();

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step);


%% Setup the initial state of the model
% The initial state of the model is setup using the model.setupInitialState() method.

initstate = model.setupInitialState();

%% Setup non-linear solver

nls = NonLinearSolver();
nls.errorOnFailure = false;

model.nonlinearTolerance = 1e-5;

%% Run simulation

model.verbose = true;
[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

return

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

