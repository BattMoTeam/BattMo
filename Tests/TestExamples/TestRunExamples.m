classdef TestRunExamples < matlab.unittest.TestCase

    properties (TestParameter)

        filename = {'runBattery1D', ...
                    'runBattery2D', ...
                    'runBattery3D', ...
                    'runChen2020', ...
                    'runCR2016'};

    end

    methods

    end

    methods (Test)

        function testRunExample(test, filename)

            % FIXME Disable plotting
            set(0, 'defaultFigureVisible', 'off');
            run(filename);
            close all
            set(0, 'defaultFigureVisible', 'on');

        end

    end

end
