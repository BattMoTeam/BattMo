%% Silicon Swelling Example 
%

%% Material properties
% We load the json structure for the material properties
jsonfilename = fullfile('ParameterData'        , ...
                        'BatteryCellParameters', ...
                        'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_nmc_silicon.json');
jsonstruct_material = parseBattmoJson(jsonfilename);

%% Geometry
% We load the json structure for the geometrical properties

jsonfilename = fullfile('Examples', 'JsonDataFiles', 'geometry1d.json');
jsonstruct_geometry = parseBattmoJson(jsonfilename);

%% Control
% We load the json structure for the geometrical properties

jsonfilename = fullfile('Examples', 'JsonDataFiles', 'cc_discharge_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);

jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                               jsonstruct_geometry, ...
                               jsonstruct_control});

%%
% The structure created in the jsonstruct follows the same hierarchy as the
% fields in the JSON input file. These can be referenced by name in the
% jsonstruct. To make life easier for ourselves we define some shorthand
% names for various parts of the structure.

ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
co      = 'Coating';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
am      = 'ActiveMaterial';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
cc      = 'CurrentCollector';

output = runBatteryJson(jsonstruct, 'runSimulation', true);

model  = output.model;
states = output.states;


%% Plotting the results
% To get the results we use the matlab cellfun function to extract the
% values Control.E, Control.I and time from each timestep (cell in the cell
% array) in states. We can then plot the vectors.

ind = cellfun(@(x) ~isempty(x), states);
states = states(ind);

E = cellfun(@(x) x.Control.E, states); 
I = cellfun(@(x) x.Control.I, states);

T = cellfun(@(x) x.time, states); 

%% Plot E as a function of the time

figure()
tiledlayout('flow');

nexttile
plot(T/hour, E)
xlabel('time [hours]')
ylabel('Cell Voltage [V]')

nexttile
plot(T/hour, I)
xlabel('time [hours]')
ylabel('Cell Current [A]')


%% Plot the overpotential of the negative electrode as a function of time

for istate = 1 : numel(states)
    states{istate} = model.addVariables(states{istate});
end

%%

figure
hold on

negativeElectrodeSize = model.(ne).grid.cells.num;
L = {};

for i = 1 : negativeElectrodeSize

    eta = cellfun(@(state) state.(ne).(co).(am).(itf).eta(i), states);
    plot(T/hour, eta);
    L{end + 1} = "x = " + int2str(i);
    
end

xlabel('time / hours')
ylabel('\eta / V')
title('Overpotential in Negative Electrode')
legend(L);

%% plot of average concentration

figure
hold on

L = {};

for i = 1 : negativeElectrodeSize

    cAver = cellfun(@(state) state.(ne).(co).(am).(sd).cAverage(i), states);
    plot(T/hour, cAver);
    L{end + 1} = "x = " + int2str(i);
    
end

xlabel('time / hours')
ylabel('c / mol/m^3')
title('Average particle concentration')
legend(L);

%% plot of Porosity

figure
hold on

negativeElectrodeSize = model.(ne).grid.cells.num;
L = {};

for i = 1 : negativeElectrodeSize

    poro = cellfun(@(state) 1 - state.(ne).(co).volumeFraction(i), states);
    plot(T/hour, poro);
    L{end + 1} = "x = " + int2str(i);
    
end

xlabel('time / hour')
ylabel('porosity / 1')
title('Negative Electrode Porosity')
legend(L);

%% Plot the radius evolution

for istate = 1 : numel(states)
    states{istate} = model.evalVarName(states{istate}, {ne, co, am, sd, 'radius'});
end
    
figure
hold on

negativeElectrodeSize = model.(ne).grid.cells.num;
L = {};

for i = 1 : negativeElectrodeSize

    radius = cellfun(@(state) state.(ne).(co).(am).(sd).radius(i), states);
    plot(T/hour, radius);
    L{end + 1} = "x = " + int2str(i);
    
end

xlabel('time [hours]')
ylabel('radius / m')
title('Silicon particle radius')
legend(L);

%% plot volume fraction
    
figure
hold on

negativeElectrodeSize = model.(ne).grid.cells.num;
L = {};

for i = 1 : negativeElectrodeSize

    vf = cellfun(@(state) state.(ne).(co).volumeFraction(i), states);
    plot(T/hour, vf);
    L{end + 1} = "x = " + int2str(i);
    
end

xlabel('time [hours]')
ylabel('volume fraction')
title('volume fraction')
legend(L);

%% Total Lithium content

figure

m = [];
vols = model.(ne).(co).grid.cells.volumes;

for istate = 1 : numel(states)
    cMaxTot = model.(ne).(co).maximumTotalConcentration;

    state = states{istate};
    x     = state.(ne).(co).(am).(sd).x;
    
    m(end + 1) = sum(cMaxTot.*x.*vols);
    
end

F = PhysicalConstants.F;
plot(T/hour, m*F/hour);
xlabel('time [hours]')
ylabel('amount / Ah');
title('Lithium content in negative electrode')

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












