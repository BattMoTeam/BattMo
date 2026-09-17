close all

jsonstruct      = parseBattmoJson('Examples/LNMO/SiGrLnmo.json');
jsonstruct_cccv = parseBattmoJson('Examples/JsonDataFiles/cccv_control.json');

jsonstruct = removeStructFields(jsonstruct, ...
                                    {'Control', 'DRate'}        , ...
                                    {'Control', 'controlPolicy'}, ...
                                    {'Control', 'dEdtLimit'}    , ...
                                    {'Control', 'dIdtLimit'}    , ...
                                    {'Control', 'rampupTime'}   , ...
                                    {'Control', 'useCVswitch'});

jsonstruct = mergeStructs({jsonstruct, jsonstruct_cccv});

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

output = runBattery(jsonstruct);

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

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
