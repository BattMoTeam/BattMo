
clear all
close all

mrstModule add ad-core mpfa

jsonstruct = parseBattmoJson(fullfile('ParameterData', 'ParameterSets', 'Chen2020', 'chen2020_lithium_ion_battery.json'));

paramobj = BatteryInputParams(jsonstruct);

% Some shorthands used for the sub-models
ne    = 'NegativeElectrode';
pe    = 'PositiveElectrode';
am    = 'ActiveMaterial';
sd    = 'SolidDiffusion';
itf   = 'Interface';
elyte = 'Electrolyte';
sei   = 'SolidElectrodeInterface';
sr    = 'SideReaction';
            
% paramobj.(ne).(am).diffusionModelType = 'full';
% paramobj.(pe).(am).diffusionModelType = 'full';
            


gen = BatteryGenerator1D();
paramobj = gen.updateBatteryInputParams(paramobj);
model = Battery(paramobj);

cgt = ComputationalGraphTool(model);

cgt.includeNodeNames = 'Negati.*Int.*R';

[g, edgelabels] = cgt.getComputationalGraph();

figure
% h = plot(g, 'edgelabel', edgelabels, 'nodefontsize', 10);
h = plot(g, 'nodefontsize', 10);

return




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
