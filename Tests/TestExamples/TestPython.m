classdef TestPython < matlab.unittest.TestCase

    properties(TestParameter)

        module = {'checkLint',  ...
                  'validateJsonScript'};

    end

    methods(Test)

        function testSetupPythonExecutable(test)
            setupPythonExecutable();
        end

        function testSetupPythonPath(test)
            setupPythonPath();
        end

        function testLoadModule(test, module)
            loadModule(module);
        end

    end

end
