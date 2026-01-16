%% Alkaline Membrane Electrolyser

%% Setup input
% Setup the physical properties for the electrolyser using json input file :battmofile:`alkalineElectrolyser.json<Electrolyser/Parameters/alkalineElectrolyser.json>`

jsonstruct_material = parseBattmoJson('Electrolyser/Parameters/alkalineElectrolyser.json');
jsonstruct_geometry = parseBattmoJson('Electrolyser/Parameters/electrolysergeometry1d.json');

jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                               jsonstruct_geometry});

simsetup = setupElectrolyserSimulation(jsonstruct);

%%
% We use the function :battmo:`rampupControl` to increase the current linearly in time

total = 10*hour;

n   = 100;
dt  = total/n;
dts = rampupTimesteps(total, dt, 5);

controlI = -3*ampere/(centi*meter)^2; % if negative, O2 and H2 are produced
tup      = total;
srcfunc  = @(time) rampupControl(time, tup, controlI, 'rampupcase', 'linear');
control  = struct('src', srcfunc);

step = struct('val', dts, 'control', ones(numel(dts), 1));
schedule = struct('control', control, 'step', step);

simsetup.schedule = schedule;

simsetup.model.verbose = true;

%% Run the simulation

states = simsetup.run();

%% Visualize the results
%
% The results contain only the primary variables of the system (the unknwons that descrive the state of the system). We
% use the method :code:`addVariables` to add all the intermediate quantities that are computed to solve the equations
% but not stored automatically in the result.

for istate = 1 : numel(states)
    states{istate} = model.addVariables(states{istate});
end

%%
% We extract the time, voltage and current values for each time step

time = cellfun(@(state) state.time, states);
E    = cellfun(@(state) state.(oer).(ptl).E, states);
I    = cellfun(@(state) state.(oer).(ctl).I, states);

%%
% We plot the results for the voltage and current

set(0, 'defaultlinelinewidth', 3)
set(0, 'defaultaxesfontsize', 15)

figure
subplot(2, 1, 1)
plot(time/hour, E)
xlabel('time [hour]');
ylabel('voltage');
title('Polarisation curve');

subplot(2, 1, 2)
plot(time/hour, -I/(1/(centi*meter)^2));
xlabel('time [hour]');
ylabel('Current [A/cm^2]');
title('Input current')

%% pH distribution plot
%
% We consider the three domains and plot the pH in each of those. We setup the helper structures to iterate over each
% domain for the plot.

models = {model.(oer).(ptl), ...
          model.(her).(ptl), ...
          model.(inm)};

fields = {{'OxygenEvolutionElectrode', 'PorousTransportLayer', 'concentrations', 2}  , ...
          {'HydrogenEvolutionElectrode', 'PorousTransportLayer', 'concentrations', 2}, ...
          {'IonomerMembrane', 'cOH'}};

h = figure();
set(h, 'position', [10, 10, 800, 450]);
hold on

ntime = numel(time);
times = linspace(1, ntime, 10);
cmap  = cmocean('deep', 10);

for ifield = 1 : numel(fields)

    fd       = fields{ifield};
    submodel = models{ifield};

    x    = submodel.grid.cells.centroids;

    for itimes = 1 : numel(times);

        itime = floor(times(itimes));
        % The method :code:`getProp` is used to recover the value from the state structure
        val   = model.getProp(states{itime}, fd);
        pH    = 14 + log10(val/(mol/litre));

        % plot of pH for the current submodel.
        plot(x/(milli*meter), pH, 'color', cmap(itimes, :));

    end

end

xlabel('x  /  mm');
ylabel('pH');
title('pH distribution in electrolyser')

colormap(cmap)
hColorbar = colorbar;
caxis([0 3]);
hTitle = get(hColorbar, 'Title');
set(hTitle, 'string', 'J (A/cm^2)');

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
