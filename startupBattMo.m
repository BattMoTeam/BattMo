% This startup file set up the MATLAB path
%
%% We use first `MRST <https://bitbucket.org/mrst/mrst-core/wiki/Home>`_ setup for MRST modules. 
% The source code for MRST is synchronized to BattMo using git-submodule mechanisms (In the MRST directory in BattMo, you
% should find the subdirectories given by the ``names`` cell array below)
% 
rootdirname = fileparts(mfilename('fullpath'));

run(fullfile(rootdirname, 'MRST', 'mrst-core',  'startup'));

names = {'autodiff', ...
         'solvers', ...
         'visualization', ...
         'model-io', ...
         'solvers'};

names = cellfun(@(x) fullfile(ROOTDIR, '..', ['mrst-', x]), names, ...
                    'UniformOutput', false);

mrstPath('addroot', names{:});

%% Register external modules

% The open source code of the 2012 version of AGMG is also available as a submodule in the directory ``Externals/agmg/``
mrstPath('register', 'agmg', fullfile(rootdirname, 'Externals', 'agmg'));

% The UPR module, which is used for meshing, is available in ``Externals/upr``
mrstPath('register', 'upr', fullfile(rootdirname, 'Externals', 'upr'));

%% The BattMo source code directories are now added directly to path

dirnames = {'Battery', 'Electrochemistry', 'Examples', 'ParameterData', 'Physics', 'Utilities', 'Tests'};

for ind = 1 : numel(dirnames)
    dirname = fullfile(rootdirname, dirnames{ind});
    addpath(genpath(dirname));
end

%% Octave requires some extra functionality
if mrstPlatform('octave')

    % Octave MRST settings
    run(fullfile('MRST', 'mrst-core', 'utils', 'octave_only', 'startup_octave.m'));

    % Disable warnings
    warning('off', 'Octave:possible-matlab-short-circuit-operator');
    warning('off', 'Octave:data-file-in-path');
    
    % Install package for json files for older octave
    if compare_versions(version, "6.4", "<=")
        try
            pkg load jsonstuff
        catch
            pkg install https://github.com/apjanke/octave-jsonstuff/releases/download/v0.3.3/jsonstuff-0.3.3.tar.gz
            pkg load jsonstuff
        end
    end

end


%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
