import matlab.unittest.TestSuite;
import matlab.unittest.selectors.HasParameter;
import matlab.unittest.parameters.Parameter;

doCreateReferenceData = false;
if doCreateReferenceData
    param = Parameter.fromData("createReferenceData", {true});
    suite = TestSuite.fromClass(?TestBattery1D, 'ExternalParameters', param);
else
    suite = TestSuite.fromClass(?TestBattery1D);
end

suite = suite.selectIf(HasParameter('Property', 'testSize', 'Value', 'short'));

suite.run()
