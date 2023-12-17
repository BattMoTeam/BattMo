clear all
close all

import matlab.unittest.TestSuite;
import matlab.unittest.selectors.HasParameter;
import matlab.unittest.parameters.Parameter;
import matlab.unittest.TestRunner

doCreateReferenceData = false;
if doCreateReferenceData
    % we create a test suite with the parameter 'createReferenceData' that is set to true, as an external parameter
    param = Parameter.fromData("createReferenceData", {true});
    suite = TestSuite.fromClass(?TestBatteryP2D, 'ExternalParameters', param);
    % This suite can then be run
    suite.run()
end

doRunTestInParallel = false;
if doRunTestInParallel
    % We run the test suite in parallel
    param = Parameter.fromData('compareWithReferenceData', {false});
    suite = TestSuite.fromClass(?TestBatteryP2D, 'ExternalParameters', param);
    suite = suite.selectIf(HasParameter('Property', 'controlPolicy', 'Value', 'CCDischarge'));
    suite = suite.selectIf(HasParameter('Property', 'use_thermal', 'Value', false));
    suite = suite.selectIf(HasParameter('Property', 'include_current_collectors', 'Value', false));
    runner = TestRunner.withTextOutput();
    result = runner.run(suite);
end

doRunFilteredSuite = false;
if doRunFilteredSuite
    % We using a subset of a suite based on parameter selection
    suite = TestSuite.fromClass(?TestBatteryP2D);
    suite = suite.selectIf(HasParameter('Property', 'testSize', 'Value', 'short'));
    suite.run();
end

doStopOnFirstFailure = false;
if doStopOnFirstFailure
   import matlab.unittest.plugins.StopOnFailuresPlugin
   suite = testsuite('TestBatteryP2D');
   runner = testrunner('textoutput');
   runner.addPlugin(StopOnFailuresPlugin)
   results = runner.run(suite);
end

doTestExamples = true;
if doTestExamples
    suite = testsuite('TestRunExamples');
    runner = testrunner('textoutput');

    import matlab.unittest.plugins.StopOnFailuresPlugin
    runner.addPlugin(StopOnFailuresPlugin)
    results = runner.run(suite);

    %results = runner.runInParallel(suite);
end


%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
