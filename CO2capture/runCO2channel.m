
%% prepare input for co2membrane
%

filename = 'CO2capture/CO2capture.json';
jsonstruct = parseBattmoJson(filename);

inputparams = CO2captureInputParams(jsonstruct);

gen = CO2captureGridGenerator();

inputparams = gen.updateInputParams(inputparams);

%% Setup model for membrane channel (using Feed input data)
%

model = CO2captureChannel(inputparams.Feed);
model.isRootSimulationModel = true;

shortNames =  {'Boundary' ,'bd'  ;
               'bcMolFractionDefinitions', 'bcmfdef';
               'bcFluxDefinition', 'bcfluxdef'};

model = model.equipModelForComputation('shortNames', shortNames);

%%
% plot computational graph
%

doplotgraph = false;

if doplotgraph
    
    cgit = model.cgit;

    close all
    cgit.plot

end

%%
% setup simulation
%

model.verbose = true;

schedule = model.Control.setupSchedule(jsonstruct);

simInput = struct('model', model, ...
                  'schedule', schedule);

simsetup = SimulationSetup(simInput);

return

[~, states] = simulateScheduleAD(initstate, model, schedule);

state = states{end};
