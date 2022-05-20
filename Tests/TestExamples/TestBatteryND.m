classdef TestBatteryND < matlab.unittest.TestCase

    properties (TestParameter)

        jsonfile = {'ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json'};
        dim = {2, 3};
    end

    methods

        function test = TestBatteryND()
            mrstModule reset
            mrstModule add ad-core mrst-gui mpfa
        end

        function states = testnd(test, dim, jsonfile)

            jsonstruct = parseBattmoJson(jsonfile);
            paramobj = BatteryInputParams(jsonstruct);

            switch dim
              case 2
                gen = BatteryGenerator2D();
                paramobj = gen.updateBatteryInputParams(paramobj);
                paramobj.NegativeElectrode.CurrentCollector.EffectiveElectricalConductivity = 1e5;
                paramobj.PositiveElectrode.CurrentCollector.EffectiveElectricalConductivity = 1e5;

              case 3
                gen = BatteryGenerator3D();
                paramobj = gen.updateBatteryInputParams(paramobj);

              otherwise
                error('dim should be only 2 or 3, not %d\n', dim);
            end

            model = Battery(paramobj);

            C      = computeCellCapacity(model);
            CRate  = 1;
            inputI = (C/hour)*CRate;

            n         = 25;
            dt        = [];
            dt        = [dt; repmat(0.5e-4, n, 1).*1.5.^[1 : n]'];
            totalTime = 1.4*hour/CRate;
            n         = 40;
            dt        = [dt; repmat(totalTime/n, n, 1)];
            times     = [0; cumsum(dt)];
            tt        = times(2 : end);
            step      = struct('val', diff(times), 'control', ones(numel(tt), 1));

            tup = 0.1;
            srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                        model.Control.Imax, ...
                                                        model.Control.lowerCutoffVoltage);
            control = struct('src', srcfunc, 'IEswitch', true);
            schedule = struct('control', control, 'step', step);

            initstate = model.setupInitialState();

            nls = NonLinearSolver();
            nls.maxIterations = 10;
            nls.errorOnFailure = true;
            nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', {{'Control', 'E'}}, ...
                                                               'targetChangeAbs', 0.03);
            model.nonlinearTolerance = 1e-5;
            model.verbose = false;

            [~, states] = simulateScheduleAD(initstate, model, schedule, ...
                                             'OutputMinisteps', true, ...
                                             'NonLinearSolver', nls);

        end

    end

    methods (Test)

        function testBattery(test, dim, jsonfile)

            states = testnd(test, dim, jsonfile);

            verifyStruct(test, states{end}, sprintf('TestBattery%dD', dim));

        end

    end

end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
