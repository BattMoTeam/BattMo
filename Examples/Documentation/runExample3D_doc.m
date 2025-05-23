%% BattMo example Json input
% This script shows an example where we setup a simulation using exclusively json input files.

%% Setting up the environment
% BattMo uses functionality from :mod:`MRST <MRSTBattMo>`. This functionality
% is collected into modules where each module contains code for doing
% specific things. To use this functionality we must add these modules to
% the matlab path by running:

mrstModule add ad-core mrst-gui

%% We load the json files
% When loading a json file using :code:`parseBattmoJson`, the output is the standard matlab structure that is
% obtained by the native matlab command :code:`jsondecode`, see `here <https://se.mathworks.com/help/matlab/ref/jsondecode.html>`_

%% Material properties
% We load the json structure for the material properties
jsonfilename = fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_nmc_graphite.json');
jsonstruct_material = parseBattmoJson(jsonfilename);

%% Geometry
% We load the json structure for the geometrical properties
jsonfilename = fullfile('Examples', 'JsonDataFiles', 'geometry3d.json');
jsonstruct_geometry = parseBattmoJson(jsonfilename);

%% Control
% We load the json structure for the geometrical properties
jsonfilename = fullfile('Examples', 'JsonDataFiles', 'cc_discharge_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);

%% Simulation parameters
% We load the json structure for the simulation parameters
jsonfilename = fullfile('Examples', 'JsonDataFiles', 'simulation_parameters.json');
jsonstruct_simparams = parseBattmoJson(jsonfilename);

%% Ouput specificiations
% We load the json structure for output extra specifications.
jsonfilename = fullfile('Examples', 'JsonDataFiles', 'extra_output.json');
jsonstruct_output = parseBattmoJson(jsonfilename);

%%
% We merge the json structures. The function issues a warning if a parameter is set with different values in the given
% structures. The rule is that the first value takes precedence.
jsonstruct = mergeJsonStructs({jsonstruct_geometry , ...
                               jsonstruct_material , ...
                               jsonstruct_control  , ...
                               jsonstruct_simparams, ...
                               jsonstruct_output   , ...
                              });


%%
% We adjust the total time with respect to the given DRate.

DRate = jsonstruct.Control.DRate;
jsonstruct.TimeStepping.totalTime = 1.4*hour/DRate;
jsonstruct.TimeStepping.numberOfTimeSteps = 40;

%% We start the simulation
% We use the function :code:`runBatteryJson` to run the simulation with json input structure

dorun = true;
if dorun
    output = runBatteryJson(jsonstruct);
end

%% Plotting

states = output.states;

E = cellfun(@(x) x.Control.E, states);
I = cellfun(@(x) x.Control.I, states);
time = cellfun(@(x) x.time, states);

figure()
subplot(1,2,1)
plot(time/hour, E)
xlabel('time [hours]')
ylabel('Cell Voltage [V]')

subplot(1,2,2)
plot(time/hour, I)
xlabel('time [hours]')
ylabel('Cell Current [A]')

figure

model = output.model
istate = 20;
plotCellData(model.grid, states{istate}.ThermalModel.T)

%%

figure

ind = 20;
state = states{ind},
E = states{ind}.Control.E
plotGrid(model.grid, 'facecolor', 'none', 'edgealpha', 0.1)
plotCellData(model.NegativeElectrode.Coating.grid, states{ind}.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface/(mol/litre))
plotCellData(model.PositiveElectrode.Coating.grid, states{ind}.PositiveElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface/(mol/litre))
title('Particle Surface Lithium Concentration');
hcb = colorbar
title(hcb, 'mol/L')
view([60, 15])

%%

T  = cellfun(@(state) max(state.ThermalModel.T), states);
figure
plot(time, T);


%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
