close all

mrstModule add ad-core

input = struct();

% If the simulation has already be run, we can just fetch the results from disk by setting runSimulation to false
runSimulation   = true;
clearSimulation = true;

input.precipitation = false;

output = run_magnesium_1D_battery(input, ...
                                  'runSimulation', runSimulation, ...
                                  'clearSimulation', clearSimulation);

model  = output.model;
states = output.states;

%% plot voltage curve

time = cellfun(@(x) x.time, states);
Enew = cellfun(@(x) x.Cathode.E, states);
% CatEnernst = cellfun(@(x) x.CathodeActiveMaterial.ENernst, states);

figure
plot((time/hour), Enew, '*-', 'linewidth', 3)
title('Potential (E)')
xlabel('time (hours)')

%% plot distributions for last time step

set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);
set(0, 'defaultlinelinewidth', 5);

simoutput = output;

casenames = {'concentrations'          , ...
             'potential'               , ...
             'quasiParticles'          , ...
             'elchemRates'};

if model.include_precipitation

    casenames = {casenames{:}       , ...
                 'precipitationRate', ...
                 'nucleation'       , ...
                 'cSat'             , ...
                 'volumeFractions'};

end

doplot = false;

if doplot

    for icase = 1 : numel(casenames)
        h = figure();
        casename = casenames{icase};
        for istate = 1 : numel(states)
            simoutput.ind = istate;
            plot1DOutputFunc(simoutput, 'casename', casename, 'fignum', h.Number);
        end
    end

end



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
