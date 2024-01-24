%% Script for running tests using github actions

% Display matlab version
disp(version)

% Debug
%delete('../Externals/mrst/mrst-core/settings.mat')

%% Display git commit ids and other stats
testdir = pwd;
[~, res] = system('git rev-parse --short HEAD');
fprintf('%s %s', pwd, res);

dirs = {'autodiff', 'core', 'model-io', 'solvers', 'visualization'};
for k = 1:numel(dirs)
    cd(sprintf('../Externals/mrst/mrst-%s', dirs{k}));
    [~, res] = system('git rev-parse --short HEAD');
    fprintf('%s %s', pwd, res);
    cd(testdir)
end

disp(pyenv)
try
    % Requires working py env, which may have not been set up
    disp(py.sys.path)
end

%% Setup BattMo
global MRST_BATCH
MRST_BATCH = true;

run('../startupBattMo.m')

%% Run python setup tests
mrstVerbose 'off';
doTestPython = false;
stopOnError = false;

if doTestPython
    suite = testsuite('TestPython');
    runner = testrunner('textoutput');

    if stopOnError
        import matlab.unittest.plugins.StopOnFailuresPlugin
        runner.addPlugin(StopOnFailuresPlugin)
    end

    results = runner.run(suite);

    % Display results
    t = table(results);
    disp(t)

    assertSuccess(results);
end

%% Run json file testing
mrstVerbose 'off';
doTestJsonFiles = true;
stopOnError = false;

if doTestJsonFiles
    suite = testsuite('TestJsonFiles');
    runner = testrunner('textoutput');

    if stopOnError
        import matlab.unittest.plugins.StopOnFailuresPlugin
        runner.addPlugin(StopOnFailuresPlugin)
    end

    results = runner.run(suite);

    % Display results
    t = table(results);
    disp(t)

    assertSuccess(results);
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
