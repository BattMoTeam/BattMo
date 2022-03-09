% This startup file set up the MATLAB path
%
%% We use first `MRST <https://bitbucket.org/mrst/mrst-core/wiki/Home>`_ setup for MRST modules. 
% The source code for MRST is synchronized to BatMo using git-submodule mechanisms (In the MRST directory in BatMo, you
% should find the subdirectories given by the ``names`` cell array below)
% 
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

%% The open source code of the 2012 version of AGMG is also available as a submodule in the directory ``Externals/agmg/``
mrstPath('register', 'agmg', fullfile(rootdirname, 'Externals/agmg/'));

%% The BatMo source code directories are now added directly to path

dirnames = {'Battery', 'Electrochemistry', 'Examples', 'ParameterData', 'Physics', 'Utilities'};

for ind = 1 : numel(dirnames)
    dirname = fullfile(rootdirname, dirnames{ind});
    addpath(genpath(dirname));
end

%% Octave requires some extra functionality
if mrstPlatform('octave')
    % Install package for json files
    try
        pkg load jsonstuff
    catch
        pkg install https://github.com/apjanke/octave-jsonstuff/releases/download/v0.3.3/jsonstuff-0.3.3.tar.gz
        pkg load jsonstuff
    end

end


%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BatMo

BatMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BatMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BatMo.  If not, see <http://www.gnu.org/licenses/>.
%}
