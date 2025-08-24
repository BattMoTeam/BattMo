%% Script for running tests using github actions
% This file must be run in the directory battmoDir()/Tests

% Display matlab version
disp(version)

% Debug
%delete('../Externals/mrst/mrst/settings.mat')

%% Display git commit ids and other stats
testdir = pwd;
[~, res] = system('git rev-parse --short HEAD');
fprintf('%s %s', pwd, res);
cd('../Externals/mrst');
[~, res] = system('git rev-parse --short HEAD');
fprintf('%s', res);
cd(testdir)

disp(pyenv)
try
    % Requires working py env, which may have not been set up
    disp(py.sys.path)
end

%% Setup BattMo
global MRST_BATCH
MRST_BATCH = true;

run('../startupBattMo.m')

% mrstSettings('set', 'allowDL', true);
% mrstSettings('set', 'promptDL', false);

import matlab.unittest.TestSuite;
import matlab.unittest.selectors.*;
import matlab.unittest.parameters.Parameter;
import matlab.unittest.TestRunner


%% Run tests

% Setup
mrstVerbose 'on';
stopOnError        = true;
runTestsInParallel = false;
doAssertSuccess    = true;

% Define which test cases to run
testCases = {
    'TestJsonFiles'  , ...
    % 'TestChen2020'   , ...
    % 'TestRunExamples', ...
    % 'TestMagnesium'
            };

% Setup test suite
for itestcase = 1 : numel(testCases)

    testCase = testCases{itestcase};
    suite =  TestSuite.fromClass(meta.class.fromName(testCase));

    if strcmp(testCase, 'TestRunExamples')

        % Tests that are not supported on github
        filenames = {'runJellyRoll'            , ...
                     'runBatteryLinearSolve'   , ...
                     'runBatteryPrecondTestP2D', ...
                     'runJuliaBridgeTest'           , ...
                     'runParameterSweep'};

        selector = HasParameter('Property', 'filename', 'Value', filenames{1});
        for ifile = 2 : numel(filenames)
            filename = filenames{ifile};
            selector = HasParameter('Property', 'filename', 'Value', filename) | selector;
        end
        suite = suite.selectIf(~selector);

    end

    suites{itestcase} = suite;

end

suite = horzcat(suites{:});

runner = testrunner('textoutput');

if stopOnError
    import matlab.unittest.plugins.StopOnFailuresPlugin;
    runner.addPlugin(StopOnFailuresPlugin);
end

% Run tests

if runTestsInParallel
    results = runner.runInParallel(suite);
else
    results = runner.run(suite);
end

t = table(results);
disp(t);

if doAssertSuccess
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
