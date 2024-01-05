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
doTestPython = true;
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
