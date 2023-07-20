function [model,schedule,state0] = setupMatlabReference(casename, jsonfolder, varargin)
%% Script for BattMo.m to produce reference solution
% - casename : is used to identify json file name for the input and saved data output
% - jsonfolder : folder where the json file is fetched
% - datafolder : folder where the computed data is saved
    
    battmo_folder = fileparts(mfilename('fullpath'));
    battmo_folder = fullfile(battmo_folder, '../..');

    % load json setup file

    json_filename = sprintf('%s.json', casename);
    json_filename = fullfile(jsonfolder, json_filename);

    if exist('parseBattmoJson') == 0
        fprintf('You need to install matlab battmo (https://github.com/BattMoTeam/BattMo).\n')
        return
    end

    jsonstruct = parseBattmoJson(json_filename);

    CRate = jsonstruct.Control.CRate;
    jsonstruct.TimeStepping.totalTime = 1.4*hour/CRate;
    jsonstruct.TimeStepping.N = 40;


    %% To run the simulation, you need to install matlab battmo

    mrstModule add ad-core mrst-gui mpfa agmg linearsolvers
    output = GenerateModelJson(jsonstruct);
    
    state0   = output.initstate;
    model    = output.model;
    schedule = output.schedule;


    %% Save solution as a matlab struct that can be imported in Julia

    if exist('class2data') == 0
        adddir = fullfile(battmo_folder, 'src/utils');
        addpath(adddir);
        fprintf('Added %s to Matlab path in order to run class2data\n', adddir);
    end

    model    = class2data(model);
    schedule = class2data(schedule);

%     filename = sprintf('%s.mat', casename);
%     filename = fullfile(datafolder, filename);
% 
%     save(filename, 'model', 'states', 'state0', "schedule")
    
%     doplot = true;
%     if doplot && runSimulation
%         ind = cellfun(@(x) not(isempty(x)), states); 
%         states = states(ind);
%         E = cellfun(@(x) x.Control.E, states); 
%         I = cellfun(@(x) x.Control.I, states);
%         time = cellfun(@(x) x.time, states);
%         plot(time, E)
%     end
end
