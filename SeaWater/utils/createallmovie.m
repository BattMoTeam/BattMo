if ~exist('dosetup')
    dosetup = true;
end

if dosetup
    close all
    states = problem.OutputHandlers.states(:);
    dosetup = false;
end

casenames = {'concentrations'   , ...
             'potential'        , ...
             'quasiparticles'   , ...
             'precipitationrate', ...
             'volumefractions'};
for icase = 1 : numel(casenames)
    casename = casenames{icase};
    createmovie(casename, model, states, 'dosavemovie', true);
end


