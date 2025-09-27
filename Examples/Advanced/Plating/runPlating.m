%% Run stand-alone active material model with lithium plating
% This example shows how to simulate a single particle of a silicon graphite
% electrode, taking into account the plating phenomenon

clear all
close all

%% Setup the properties of Li-ion battery materials and cell design

filename = fullfile('ParameterData'        , ...
                    'BatteryCellParameters', ...
                    'LithiumIonBatteryCell', ...
                    'lithium_ion_battery_nmc_graphite.json');

jsonstruct = parseBattmoJson(filename);

%%
% We define some shortcuts

ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
lp      = 'LithiumPlating';
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

jsonstruct.(ne).(co).(am).diffusionModelType = 'full';
jsonstruct.(pe).(co).(am).diffusionModelType = 'full';

jsonstruct.(ne).(co).(am).useLithiumPlating = true;

%%
% OCP is computed via a function described in the article (S-5)
jsonstruct.(ne).(co).(am).(itf).openCircuitPotential.functionName = 'computeOCP_Graphite_Latz';
jsonstruct.(ne).(co).(am).(itf).includeEntropyChange = false;

%%
% Flag for stand-alone model
%

jsonstruct.(ne).(co).(am).isRootSimulationModel = true;

jsonstruct_lithium_plating = parseBattmoJson(fullfile('Examples', 'Advanced', 'Plating', 'lithium_plating.json'));

jsonstruct.(ne).(co).(am).LithiumPlating = jsonstruct_lithium_plating.LithiumPlating;

%%
% Setup InputParams

inputparams = BatteryInputParams(jsonstruct);
inputparams = inputparams.(ne).(co).(am);


%% Setup the model

model = ActiveMaterial(inputparams);

%%
% We equip the model for simulation
%

model = model.setupForSimulation();
model.verbose = true;


%% setup schedule
% This part is essential to see the lithium plating effect. Increasing Iref strengthens the effect.

Iref = 7e-13;
Imax = Iref;
total = 6e-3*hour*(Iref/Imax); % total time of the charge
n     = 400;
dt    = total/n;
step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

tup = 1*second*(Iref/Imax);

srcfunc = @(time) rampupControl(time, tup, -Imax);

control.src = srcfunc;

schedule = struct('control', control, 'step', step);

%% Setup initial state

sd  = 'SolidDiffusion';
itf = 'Interface';

cElectrolyte   = 5e-1*mol/litre;
phiElectrolyte = 0;
T              = 298;
cElectrodeInit = 30*mol/litre;

[model, initstate] = setupPlatingInitialState(model, T, cElectrolyte, phiElectrolyte, cElectrodeInit, Imax);

%% Run first simulation, particle charge
%

inputSim = struct('model'    , model   , ...
                  'schedule' , schedule, ...
                  'initstate', initstate);
simsetup = SimulationSetup(inputSim);

states = simsetup.run();

chargeStates = states; % for later

%% Setup discharge simulation
% 
%
% We use the last state of the previous simulation to initialise the new one

simsetup.initstate = states{end};

%%
% Setup schedule structure for discharge
%

cmin = (model.(itf).guestStoichiometry0)*(model.(itf).saturationConcentration);

control.stopFunction = @(model, state, state0_inner) (state.(sd).cSurface <= cmin);
control.src          = @(time) rampupControl(time, tup, Imax);

simsetup.schedule = struct('control', control, 'step', step);

%% Run second simulation, particle discharge
%

states = simsetup.run();

%%
% We concatenate the two phases

dischargeStates = states;
states = vertcat(chargeStates, dischargeStates);

%% Plotting

ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

time     = cellfun(@(state) state.time, states);
cSurface = cellfun(@(state) state.(sd).cSurface, states);
E        = cellfun(@(state) state.E, states);

figure
plot(time, cSurface/(1/litre));
xlabel('time [second]');
ylabel('Surface concentration [mol/L]');
title('Surface concentration');

figure
plot(time, E);
xlabel('time [second]');
ylabel('Potential [mol/L]');
title('Potential difference');

cmin = cellfun(@(state) min(state.(sd).c), states);
cmax = cellfun(@(state) max(state.(sd).c), states);

for istate = 1 : numel(states)
    states{istate} = model.evalVarName(states{istate}, {sd, 'cAverage'});
end

lp = 'LithiumPlating';

varsToEval = {{'Interface'     , 'eta'}         , ...
              {'LithiumPlating', 'etaPlating'}  , ...
              {'LithiumPlating', 'etaChemical'} , ...
              {'Interface'     , 'intercalationFlux'}           , ...
              {'LithiumPlating', 'platingFlux'} , ...
              {'LithiumPlating', 'chemicalFlux'}, ...
              {'LithiumPlating', 'surfaceCoverage'}, ...
              {'LithiumPlating', 'platedThickness'}};
for k = 1:numel(states)
    for var = 1:numel(varsToEval)
        states{k} = model.evalVarName(states{k}, varsToEval{var});
    end
end

varnames = {
            'platedConcentration', ...    
            'eta', ...             
            'etaPlating', ...      
            'etaChemical', ...     
            'platingFlux', ...     
            'chemicalFlux', ...    
            'intercalationFlux', ...               
            'surfaceCoverage', ...
            'platedThickness'};

vars = {};

vsa = model.LithiumPlating.volumetricSurfaceArea;

vars{end + 1} = cellfun(@(s) s.(lp).platedConcentration, states);
vars{end + 1} = cellfun(@(s) s.(itf).eta, states);
vars{end + 1} = cellfun(@(s) s.(lp).etaPlating, states);
vars{end + 1} = cellfun(@(s) s.(lp).etaChemical, states);
vars{end + 1} = cellfun(@(s) s.(lp).platingFlux .* s.(lp).surfaceCoverage .* vsa, states);
vars{end + 1} = cellfun(@(s) s.(lp).chemicalFlux .* s.(lp).surfaceCoverage, states);
vars{end + 1} = cellfun(@(s) s.(itf).intercalationFlux .* (1 - s.(lp).surfaceCoverage), states);
vars{end + 1} = cellfun(@(s) s.(lp).surfaceCoverage, states);
vars{end + 1} = cellfun(@(s) s.(lp).platedThickness, states);

% Variable : surfaceCoverage
figure
plot(time, vars{8}, '-');
xlabel('time [second]');
ylabel(varnames{8});
title(varnames{8});

%%
% Area fraction that is plated. The more lithium is plated, the less it can
% pass from the electrolyte to intercalate into the electrode. Thus, if
% surfaceCoverage = 1, no more lithium can be intercalated from the solution. 
% The electrode can still be filled with plated lithium through the
% chemical flux

% Variable : platingFlux .* surfaceCoverage
figure
plot(time, vars{5}, '-');
xlabel('time [second]');
ylabel(varnames{5});
title(varnames{5});

%%
% We can see clearly here the 4 differents steps of the lithium plating phenomenon.
% First, at the end of the charge, the amount of lithium being plated grows
% faster and faster as the area where plating is possible increases. 
% 
% Then, as the plated lithium covers the whole particle, the plating flux stabilises. 
%
% At the beginning of the discharge, the plated lithium begins to strip, as
% the plated layer is the only electron source available (the intercalated
% lithium has no contact with the electrolyte)

% Finally, the surfaceCoverage decreases, resulting in a slower stripping.

% Variable : chemicalFlux .* surfaceCoverage
figure
plot(time, vars{6}, '-');
xlabel('time [second]');
ylabel(varnames{6});
title(varnames{6});

% Variable : intercalationFlux .* (1 - surfaceCoverage)
figure
plot(time, vars{7}, '-');
xlabel('time [second]');
ylabel(varnames{7});
title(varnames{7});
% No more lithium is intercalated during the time the whole surface is covered with plated lithium.

% Variable : platedThickness
figure
plot(time, vars{9}, '-');
xlabel('time [second]');
ylabel(varnames{9});
title(varnames{9});

figure
plot(time, vars{2}, '-');
xlabel('time [second]');
ylabel(varnames{2});
title(varnames{2});

figure
plot(time, vars{3}, '-');
xlabel('time [second]');
ylabel(varnames{3});
title(varnames{3});

%% Sum up

figure;

yyaxis left
plot(time, cSurface/(1/litre), '-', 'LineWidth', 1.5);
ylabel('Surface concentration [mol/L]');

yyaxis right
plot(time, vars{5}, '-', 'LineWidth', 1.5);
ylabel('Volumetric plating flux [mol/(m³·s)]');

xlabel('Time [second]');
title('Surface Concentration and Volumetric plating flux');
legend('Surface concentration', 'Plating flux');
grid on;

