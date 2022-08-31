
clear all
close all

mrstModule add ad-core mpfa

jsonstruct = parseBattmoJson('ParameterData/ParameterSets/Chen2020/chen2020_lithium_ion_battery.json');

paramobj = BatteryInputParams(jsonstruct);

% Some shorthands used for the sub-models
ne    = 'NegativeElectrode';
pe    = 'PositiveElectrode';
am    = 'ActiveMaterial';
sd    = 'SolidDiffusion';
itf   = 'Interface';
elyte = 'Electrolyte';

paramobj.(ne).(am).(sd).diffusionModelType = 'full';
paramobj.(pe).(am).(sd).diffusionModelType = 'full';


gen = BareBatteryGenerator3D();
paramobj = gen.updateBatteryInputParams(paramobj);
model = Battery(paramobj);

submodel = model.NegativeElectrode.ActiveMaterial.Interface.registerVarAndPropfuncNames();
% submodel = model.NegativeElectrode.ActiveMaterial.registerVarAndPropfuncNames();
% submodel = model.Electrolyte.registerVarAndPropfuncNames();

[g, edgelabels] = setupGraph(submodel);

cgf = ComputationalGraphFilter(g);
cgf.includeNodeNames = [];

g = cgf.setupGraph();

figure
h = plot(g, 'nodefontsize', 18);

doaddlabel = false;
if doaddlabel
    labeledge(h, (1 : g.numedges), edgelabels);
end





%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
