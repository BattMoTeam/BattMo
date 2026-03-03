close all

jsonstruct      = parseBattmoJson('Examples/LNMO/SiGrLnmo.json');
jsonstruct_cccv = parseBattmoJson('Examples/JsonDataFiles/cccv_control.json');

jsonstruct = removeJsonStructFields(jsonstruct, ...
                                    {'Control', 'DRate'}        , ...
                                    {'Control', 'controlPolicy'}, ...
                                    {'Control', 'dEdtLimit'}    , ...
                                    {'Control', 'dIdtLimit'}    , ...
                                    {'Control', 'rampupTime'}   , ...
                                    {'Control', 'useCVswitch'});

jsonstruct = mergeJsonStructs({jsonstruct, jsonstruct_cccv});

jsonstruct.TimeStepping.numberOfTimeSteps = 100;

jsonstruct.Control.DRate = 1;
jsonstruct.Control.CRate = 1;
jsonstruct.Control.lowerCutoffVoltage = 3.8;
jsonstruct.Control.upperCutoffVoltage = 4.8;

jsonstruct.Control.initialControl = 'charging';

jsonstruct.SOC = 0;
jsonstruct.NegativeElectrode.Coating.ActiveMaterial2.massFraction = 0.04;
jsonstruct.NegativeElectrode.Coating.ActiveMaterial1.massFraction = 0.96;
jsonstruct.Control.numberOfCycles = 1;

output = runBatteryJson(jsonstruct);

states = output.states;
model  = output.model;

%% plot


set(0, 'defaultlinelinewidth', 3);
set(0, 'defaultaxesfontsize', 15);

t = cellfun(@(state) state.time, states);
E = cellfun(@(state) state.Control.E, states);

figure
plot(t, E)

%% Shortcuts
% We define shorcuts for the sub-models.

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


figure
hold on

for istate = 1 : numel(states)
    states{istate} = model.evalVarName(states{istate}, {ne, co, 'SOC'});
end

SOC  = cellfun(@(x) x.(ne).(co).SOC, states);
SOC1 = cellfun(@(x) x.(ne).(co).(am1).SOC, states);
SOC2 = cellfun(@(x) x.(ne).(co).(am2).SOC, states);

plot(t/hour, SOC, 'displayname', 'SOC - cumulated');
plot(t/hour, SOC1, 'displayname', 'SOC - Graphite');
plot(t/hour, SOC2, 'displayname', 'SOC - Silicon');

xlabel('Time / h');
ylabel('SOC / -');
title('SOCs')

legend show

