
rootdirname = fileparts(mfilename('fullpath'));

run(fullfile(rootdirname, 'MRST/mrst-core/startup'));

names = {'autodiff', ...
         'solvers', ...
         'visualization', ...
         'model-io', ...
         'solvers'};

names = cellfun(@(x) fullfile(ROOTDIR, '..', ['mrst-', x]), names, ...
                    'UniformOutput', false);

mrstPath('addroot', names{:});

mrstPath('register', 'multimodel', fullfile(rootdirname, 'MRST/mrst-multimodel/'));

dirnames = {'Battery', 'Electrochemistry', 'Examples', 'Materials', 'Physics', 'Utilities'};

for ind = 1 : numel(dirnames)
    dirname = fullfile(rootdirname, dirnames{ind});
    addpath(genpath(dirname));
end
