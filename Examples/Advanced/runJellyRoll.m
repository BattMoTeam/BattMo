clear
close all

%% We load the geometrical parameters for a jellyroll (4680 model)

jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', '4680-geometry.json'));

testing = true;
if testing
    % We setup a smaller model for quicker testing
    fprintf('We setup a smaller case for quicker testing\n');
    rOuter = jsonstruct_geometry.Geometry.innerRadius + 1*milli*meter;
    jsonstruct_geometry.Geometry.outerRadius                         = rOuter;
    jsonstruct_geometry.Geometry.numberOfDiscretizationCellsVertical =  2;

    tabparams = struct('usetab', false);
    jsonstruct_geometry.NegativeElectrode.CurrentCollector.tabparams = tabparams;
    jsonstruct_geometry.PositiveElectrode.CurrentCollector.tabparams = tabparams;
    
end


%% We load material parameters

jsonstruct_material = parseBattmoJson(fullfile('ParameterData'        , ...
                                               'BatteryCellParameters', ...
                                               'LithiumIonBatteryCell', ...
                                               'lithium_ion_battery_nmc_graphite.json'));

jsonstruct_material = removeJsonStructField(jsonstruct_material, {'include_current_collectors'});

%% We load the control parameters

jsonstruct_control = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'cc_discharge_control.json'));

%% We merge all the parameters

jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                               jsonstruct_geometry, ...
                               jsonstruct_control});

jsonstruct.include_current_collectors = true;

output = runBatteryJson(jsonstruct);

%% Process output and recover the output voltage and current from the output states.

states = output.states;

E = cellfun(@(x) x.Control.E, states);
I = cellfun(@(x) x.Control.I, states);
time = cellfun(@(x) x.time, states);

figure
plot(time, E, 'linewidth', 3);
set(gca, 'fontsize', 18);
title('Cell Voltage / V')
xlabel('time')


%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
