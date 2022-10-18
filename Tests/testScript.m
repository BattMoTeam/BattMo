import matlab.unittest.TestSuite;
import matlab.unittest.selectors.HasParameter;
import matlab.unittest.parameters.Parameter;
import matlab.unittest.TestRunner

doCreateReferenceData = false;
if doCreateReferenceData
    % we create a test suite with the parameter 'createReferenceData' that is set to true, as an external parameter
    param = Parameter.fromData("createReferenceData", {true});
    suite = TestSuite.fromClass(?TestBattery1D, 'ExternalParameters', param);
    % This suite can then be run 
    suite.run()
end

doRunTestInParallel = false;
if doRunTestInParallel
    % We run the test suite in parallel
    suite = TestSuite.fromClass(?TestBattery1D);
    runner = TestRunner.withNoPlugins();
    result = runner.runInParallel(suite);
end

doRunFilteredSuite = false;
if doRunFilteredSuite
    % We using a subset of a suite based on parameter selection
    suite = TestSuite.fromClass(?TestBattery1D);
    suite = suite.selectIf(HasParameter('Property', 'testSize', 'Value', 'short'));
    suite.run();
end

doStopOnFirstFailure = false;
if doStopOnFirstFailure
   import matlab.unittest.plugins.StopOnFailuresPlugin
   suite = testsuite('TestBattery1D');
   runner = testrunner('textoutput');
   runner.addPlugin(StopOnFailuresPlugin)
   results = runner.run(suite);
end
