%Startup JuliaBridge by adding important folders to path
if ~exist('setupMatlabModel')
        adddir = fullfile(battmo_folder, '..',GenerateModel');
        addpath(adddir);
        fprintf('Added %s to Matlab path in order to run setupMatlabModel', adddir);
end

if ~exist('ServerManager')
        adddir = fullfile(battmo_folder, '..','JuliaInterface');
        addpath(adddir);
        fprintf('Added %s to Matlab path in order to run ServerManager', adddir);
end