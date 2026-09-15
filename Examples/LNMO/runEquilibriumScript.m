close all
mrstDebug(20);

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


ne  = 'NegativeElectrode';
pe  = 'PositiveElectrode';
co  = 'Coating';
am1 = 'ActiveMaterial1';
am2 = 'ActiveMaterial2';
am  = 'ActiveMaterial';
itf = 'Interface';

jsonstruct.(ne).(co).(am1).(itf).guestStoichiometry100 = 1;
jsonstruct.(ne).(co).(am1).(itf).guestStoichiometry0   = 0.01;
jsonstruct.(ne).(co).(am2).(itf).guestStoichiometry100 = 1;
jsonstruct.(ne).(co).(am2).(itf).guestStoichiometry0   = 0.01;

jsonstruct.(pe).(co).(am).(itf).guestStoichiometry0 = 1;

css = CellSpecificationSummary(jsonstruct);

css.printSpecifications();

return


model = setupModelFromJson(jsonstruct);

ecs = EquilibriumConcentrationSolver(model);

[E, extras] = advancedComputeCellEnergy(model);
m = computeCellMass(model);

fprintf('Specific Energy %g\n', (E/m)/hour);

dsfunc = extras.dischargeFunction;
soc = linspace(0, 1, 100);

figure
plot(soc, dsfunc(soc));
xlabel('SOC / -');
ylabel('Voltage / V');

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
