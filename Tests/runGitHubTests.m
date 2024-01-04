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

%% Run tests

mrstVerbose 'on';

doTestJsonFiles = true;
stopOnError = true;

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




% %%
% mrstVerbose 'on';
% mrstSettings('set', 'useMEX', false);
% mrstModule add ad-core mpfa
% addpath('../MRST-debug')

% % Setup tests
% import matlab.unittest.TestSuite;
% import matlab.unittest.selectors.HasParameter;
% import matlab.unittest.parameters.Parameter;
% import matlab.unittest.TestRunner

% % Run tests
% %suite = TestSuite.fromFolder('TestExamples');
% %suite = suite.selectIf(HasParameter('Property', 'testSize', 'Value', 'short'));
% %suite = TestSuite.fromClass(?TestBattery1D);

% param = Parameter.fromData("createReferenceData", {false});
% suite = TestSuite.fromClass(?TestBattery1D, 'ExternalParameters', param);

% params = {'Property', 'controlPolicy', 'Value', 'CCCV',...
%           'Property', 'use_thermal', 'Value', true,...
%           'Property', 'include_current_collectors', 'Value', true,...
%           'Property', 'diffusionModelType', 'Value', 'simple',...
%           'Property', 'testSize', 'Value', 'short'};
% for k = 1:numel(params)/4
%     p = params((4*k-3):(4*k));
%     suite = suite.selectIf(HasParameter(p{:}));
% end

% results = suite.run();

% % Display results
% t = table(results)

% % Assert
% assertSuccess(results);
