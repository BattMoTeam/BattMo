% clear all
% close all




% We define some shorthand names for simplicity.
ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
co      = 'Coating';
am      = 'ActiveMaterial';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
cc      = 'CurrentCollector';

jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));

jsonstruct = mergeStructs({jsonstruct_material, ...
                               jsonstruct_geometry});

includeDoubleLayer = false;

if includeDoubleLayer

    jsonstruct.(ne).(co).(am).(itf).useDoubleLayerCapacity = true;
    jsonstruct.(ne).(co).(am).(itf).doubleLayerCapacitance = 0.2;

end

[model, inputparams, ~, gen] = setupModelFromJson(jsonstruct);

c_ne = 29.866*mol/litre; % initial concentration at negative electrode
c_pe = 17.038*mol/litre; % initial concentration at positive electrode

initstate = initStateChen2020(model, c_ne, c_pe);

options = [];
options.stateInitialization.initializationSetup = 'given state';
options.stateInitialization.computeSteadyState  = false;

extrastructs = [];
extrastructs.initstate = initstate;

impsolv = ImpedanceSolver(inputparams, options, extrastructs);







set(0, 'defaultlinelinewidth', 3);



figure


omegas = logspace(-4, 2, 30);
Z = impsolv.computeImpedance(omegas);
hold on

plot(real(Z), -imag(Z), 'displayname', 'battmo');

axis equal;


docompare = true;
if docompare
    p = fileparts(mfilename('fullpath'));
    data = load(fullfile(p, 'utils', 'pybamm_chen_impedances.mat'));
    Zpybamm = data.impedances;
    plot(real(Zpybamm), -imag(Zpybamm), 'displayname', 'pybamm');
end

xlabel('real(Z) / Ω')
ylabel('-imag(Z) / Ω')
legend show
title('Impedance')

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
