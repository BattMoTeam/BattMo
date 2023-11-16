function export = setupMatlabModel(casename, jsonfolder, generate_reference_solution,varargin)
%% Script for BattMo.m to produce reference solution
% - casename : is used to identify json file name for the input and saved data output
% - jsonfolder : folder where the json file is fetched
% - datafolder : folder where the computed data is saved
    
    % load json setup file

    json_filename = sprintf('%s.json', casename);
    json_filename = fullfile(jsonfolder, json_filename);

    jsonstruct = parseBattmoJson(json_filename);

    CRate = jsonstruct.Control.CRate;
    jsonstruct.TimeStepping.totalTime = 1.4*hour/CRate;
    jsonstruct.TimeStepping.numberOfTimeSteps = 100;
    
    %% To run the simulation, you need to install matlab battmo

    export = runBatteryJson(jsonstruct, 'runSimulation', generate_reference_solution);

    %% Save solution as a matlab struct that can be imported in Julia

    export.model    = class2data(export.model);
    export.schedule = class2data(export.schedule);
    
end
