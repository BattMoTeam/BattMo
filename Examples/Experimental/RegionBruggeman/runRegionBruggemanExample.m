%% Setup material properties

jsonfilename = fullfile('Examples', 'Experimental', 'RegionBruggeman', 'region_bruggeman.json');
jsonstruct_bruggeman = parseBattmoJson(jsonfilename);

jsonfilename = fullfile('ParameterData'        , ...
                        'BatteryCellParameters', ...
                        'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_lco_graphite.json');
jsonstruct_material = parseBattmoJson(jsonfilename);

jsonstruct_material.include_current_collectors = true;
jsonstruct_material.use_thermal = false;

%% Setup geometry
jsonfilename = fullfile('Examples'     , ...
                        'JsonDataFiles', ...
                        'geometry1d.json');
jsonstruct_geometry = parseBattmoJson(jsonfilename);

%% Setup Control
jsonfilename = fullfile('Examples', 'JsonDataFiles', 'cc_discharge_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);

jsonstruct_timestepping.TimeStepping.numberOfTimeSteps = 40;

jsonstruct = mergeStructs({jsonstruct_bruggeman, ...
                               jsonstruct_geometry , ...
                               jsonstruct_material , ...
                               jsonstruct_control  , ...
                               jsonstruct_timestepping});

output = runBattery(jsonstruct);

doplot = true;

if doplot

    states = output.states;
    E = cellfun(@(state) state.Control.E, states);
    t = cellfun(@(state) state.time, states);

    figure
    plot(t, E)
    
end

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
