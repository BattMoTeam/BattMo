%% Simulation of a composite active material
% This example shows how to simulate a composite active material

%% Setup the properties of the battery
%
% We load the property of a composite silicon graphite electrode. 

jsonstruct_composite_material = parseBattmoJson('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/composite_silicon_graphite.json');

%%
% In this structure, we have the material property of two active materials, |ActiveMaterial1| and |ActiveMaterial2|.
flattenJsonStruct(jsonstruct_composite_material);

%%
% For the remaining properties, we load a standard data set
jsonstruct_cell = parseBattmoJson('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json');

%%
% We remove from this structure active material field. This step is not necessary but is cleaner and we avoid a
% warning.
jsonstruct_cell = removeJsonStructField(jsonstruct_cell, {'NegativeElectrode', 'Coating', 'ActiveMaterial'});

%%
% we load a 1d geometry
jsonfilename = fullfile('Examples', 'JsonDataFiles', 'geometry1d.json');
jsonstruct_geometry = parseBattmoJson(jsonfilename);

%%
% We merge the json structures
jsonstruct = mergeJsonStructs({jsonstruct_composite_material, ...
                               jsonstruct_cell              , ...
                               jsonstruct_geometry});

%%
% We do not consider the thermal model and remove the current collector. We also use a CV switch control.
jsonstruct.use_thermal                = false;
jsonstruct.include_current_collectors = false;

%% 
% We define some shorcuts for the sub-models, for convenience

ne   = 'NegativeElectrode';
pe   = 'PositiveElectrode';
co   = 'Coating';
am1  = 'ActiveMaterial1';
am2  = 'ActiveMaterial2';
bd   = 'Binder';
ad   = 'ConductingAdditive';
sd   = 'SolidDiffusion';
itf  = 'Interface';
ctrl = 'Control';

%% We modify some parameters
% We adjust the mass fractions  parameters of the active material in the negative electrode

jsonstruct = setJsonStructField(jsonstruct, {ne, co, am1, 'massFraction'}, 0.9, 'handleMisMatch', 'quiet');
jsonstruct = setJsonStructField(jsonstruct, {ne, co, am2, 'massFraction'}, 0.08, 'handleMisMatch', 'quiet');
jsonstruct = setJsonStructField(jsonstruct, {ne, co, bd , 'massFraction'}, 0.01, 'handleMisMatch', 'quiet');
jsonstruct = setJsonStructField(jsonstruct, {ne, co, ad , 'massFraction'}, 0.01, 'handleMisMatch', 'quiet');

%% We run the simulations
%
output = runBatteryJson(jsonstruct);

%% Plotting
% We extract the voltage, current and time from the simulation output

states = output.states;
model  = output.model;

E    = cellfun(@(x) x.Control.E, states);
I    = cellfun(@(x) x.Control.I, states);
time = cellfun(@(x) x.time, states);

%%
% We plot the voltage and current

figure
subplot(2, 1, 1);
plot(time/hour, E);
xlabel('Time / h');
ylabel('Voltage / V');
title('Voltage')
subplot(2, 1, 2);
plot(time/hour, I/milli);
xlabel('Time / h');
ylabel('Current / mA');
title('Current')

%%
% We compute and plot the state of charges in the different material

figure
hold on

for istate = 1 : numel(states)
    states{istate} = model.evalVarName(states{istate}, {ne, co, 'SOC'});
end

SOC  = cellfun(@(x) x.(ne).(co).SOC, states);
SOC1 = cellfun(@(x) x.(ne).(co).(am1).SOC, states);
SOC2 = cellfun(@(x) x.(ne).(co).(am2).SOC, states);

plot(time/hour, SOC, 'displayname', 'SOC - cumulated');
plot(time/hour, SOC1, 'displayname', 'SOC - Graphite');
plot(time/hour, SOC2, 'displayname', 'SOC - Silicon');

xlabel('Time / h');
ylabel('SOC / -');
title('SOCs')

legend show

%% plot of the particle concentration distribution in the particle at the end time
%
%%
% We recover the state at the last time step
%

state = states{end};

%%
% We iterate over the two active materials. The first one is the graphite and the second one the silicon
%
ams = {am1, am2};

for iam = 1 : numel(ams)

    am = ams{iam};

    model_sd = model.(ne).(co).(am).(sd);
    state_sd = state.(ne).(co).(am).(sd);

    %%
    % We recover the concentration as an array. The column index corresponds to the spatial direction (here x as we are considering a 1D model) and the row index corresponds to the particle radia direction
    c = model_sd.getParticleConcentrations(state_sd);
    r = model_sd.operators.radii;

    figure
    hold on
    %%
    % We plot the concentration distribution at the last point in the grid, which corresponds in this case to the
    % closest to the positive electrode
    plot(r/(micro*meter), c(size(c, 1), :)/(mol/litre));
    title(sprintf('Particle concentration profile in %s', am));
    xlabel('radius / m');
    ylabel('concentration / mol/litre');
    
end


%% Charge step

initstate = states{end};
jsonstruct.(ctrl).CRate = 1;
jsonstruct = setJsonStructField(jsonstruct, {'Control', 'controlPolicy'}, 'CCCharge', 'handleMisMatch', 'quiet');

jsonstruct.initializationSetup = 'given matlab object';

output = runBatteryJson(jsonstruct, 'initstate', initstate);

%% Visualisation

%%
% We concatenate the states we have computed

dischargeStates = states; % see previous assignement
chargeStates    = output.states; % assigned from last simulation output

allStates = vertcat(dischargeStates, chargeStates);

%%

% We extract the voltage, current and time from the simulation output
E    = cellfun(@(x) x.Control.E, allStates);
I    = cellfun(@(x) x.Control.I, allStates);
time = cellfun(@(x) x.time, allStates);

%%
%  We plot the voltage and current
figure
subplot(2, 1, 1);
plot(time/hour, E);
xlabel('Time / h');
ylabel('Voltage / V');
title('Voltage')
subplot(2, 1, 2);
plot(time/hour, I/milli);
xlabel('Time / h');
ylabel('Current / mA');
title('Current')

%%
% We compute and plot the state of charges in the different material

figure
hold on

for istate = 1 : numel(allStates)
    allStates{istate} = model.evalVarName(allStates{istate}, {ne, co, 'SOC'});
end

SOC  = cellfun(@(x) x.(ne).(co).SOC, allStates);
SOC1 = cellfun(@(x) x.(ne).(co).(am1).SOC, allStates);
SOC2 = cellfun(@(x) x.(ne).(co).(am2).SOC, allStates);

plot(time/hour, SOC, 'displayname', 'SOC - cumulated');
plot(time/hour, SOC1, 'displayname', 'SOC - Graphite');
plot(time/hour, SOC2, 'displayname', 'SOC - Silicon');

xlabel('Time / h');
ylabel('SOC / -');
title('SOCs')

legend show

