close all

mrstModule add ad-core

input = struct();

% If the simulation has already be run, we can just fetch the results from disk by setting runSimulation to false
runSimulation = false;
clearSimulation = false;

input.precipitation = false;

output = run_magnesium_1D_battery(input, 'runSimulation', runSimulation, 'clearSimulation', clearSimulation);

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

for icase = 1 : numel(casenames)
    h = figure();
    casename = casenames{icase};
    for istate = 1 : numel(states)
        simoutput.ind = istate;
        plot1DOutputFunc(simoutput, 'casename', casename, 'fignum', h.Number);
    end
end


