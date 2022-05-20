classdef TestChen2020 < matlab.unittest.TestCase

    properties (TestParameter)
        
        useSimplifiedDiffusionPE = {true, false};
        useSimplifiedDiffusionNE = {true, false};

    end
    
    methods

        function test = TestChen2020()
            mrstModule reset
            mrstModule add ad-core mrst-gui mpfa
        end
        
        function states = testchen2020(test, doSingleTimestep, useSimplifiedDiffusionPE, useSimplifiedDiffusionNE)

            jsonstruct = parseBattmoJson('ParameterData/ParameterSets/Chen2020/chen2020_lithium_ion_battery.json');
            paramobj = BatteryInputParams(jsonstruct);

            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            sd    = 'SolidDiffusion';
            elyte = 'Electrolyte';
            itf   = 'Interface';
            ctrl  = 'Control';

            % Battery modeland params
            gen = BareBatteryGenerator3D();
            paramobj = gen.updateBatteryInputParams(paramobj);
            paramobj.(ne).(am).InterDiffusionCoefficient = 0;
            paramobj.(pe).(am).InterDiffusionCoefficient = 0;
            paramobj.(ne).(am).(sd).useSimplifiedDiffusionModel = useSimplifiedDiffusionNE;
            paramobj.(pe).(am).(sd).useSimplifiedDiffusionModel = useSimplifiedDiffusionPE;

            model = Battery(paramobj);
            
            model.Control.Imax = 5;

            % Schedule
            fac   = 2; 
            total = 1.4*hour; 
            n     = 100;
            dt0   = total*1e-6; 
            times = getTimeSteps(dt0, n, total, fac); 
            dt    = diff(times);
            if doSingleTimestep
                dt = dt(1);
            else
                dt = dt(1:end);
            end
            step = struct('val', dt, 'control', ones(size(dt)));

            tup = 0.1;
            srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                        model.Control.Imax, ...
                                                        model.Control.lowerCutoffVoltage);
            control = struct('src', srcfunc, 'IEswitch', true);
            schedule = struct('control', control, 'step', step); 

            % Initial state
            nc = model.G.cells.num;
            T = model.initT;
            initstate.ThermalModel.T = T*ones(nc, 1);

            bat = model;
            initstate = model.updateTemperature(initstate);

            % setup negative electrode initial state
            nitf = model.(ne).(am).(itf); 

            c = 29866.0;
            if model.(ne).(am).useSimplifiedDiffusionModel
                nenp = model.(ne).(am).G.cells.num;
                initstate.(ne).(am).c = c*ones(nenp, 1);
            else
                nenr = model.(ne).(am).(sd).N;
                nenp = model.(ne).(am).(sd).np;
                initstate.(ne).(am).(sd).c = c*ones(nenr*nenp, 1);
            end

            initstate.(ne).(am).(sd).cSurface = c*ones(nenp, 1);

            initstate.(ne).(am) = model.(ne).(am).updateConcentrations(initstate.(ne).(am));
            initstate.(ne).(am).(itf) = nitf.updateOCP(initstate.(ne).(am).(itf));

            OCP = initstate.(ne).(am).(itf).OCP;
            ref = OCP(1);

            initstate.(ne).(am).phi = OCP - ref;

            % setup positive electrode initial state
            pitf = bat.(pe).(am).(itf); 

            c = 17038.0;

            if model.(pe).(am).useSimplifiedDiffusionModel
                penp = model.(pe).(am).G.cells.num;
                initstate.(pe).(am).c = c*ones(penp, 1);
            else
                penr = model.(pe).(am).(sd).N;
                penp = model.(pe).(am).(sd).np;
                initstate.(pe).(am).(sd).c = c*ones(penr*penp, 1);
            end
            initstate.(pe).(am).(sd).cSurface = c*ones(penp, 1);

            initstate.(pe).(am) = model.(pe).(am).updateConcentrations(initstate.(pe).(am));
            initstate.(pe).(am).(itf) = pitf.updateOCP(initstate.(pe).(am).(itf));

            OCP = initstate.(pe).(am).(itf).OCP;

            initstate.(pe).(am).phi = OCP - ref;

            initstate.(elyte).phi = zeros(bat.(elyte).G.cells.num, 1) - ref;
            initstate.(elyte).c = 1000*ones(bat.(elyte).G.cells.num, 1);

            % setup initial positive electrode external coupling values
            initstate.(ctrl).E = OCP(1) - ref;
            initstate.(ctrl).I = 0;
            initstate.(ctrl).ctrlType = 'constantCurrent';

            nls = NonLinearSolver(); 
            nls.maxIterations = 10; 
            nls.errorOnFailure = true;
            model.nonlinearTolerance = 1e-5; 
            model.verbose = false;

            model.AutoDiffBackend = AutoDiffBackend();

            [~, states] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', false, 'NonLinearSolver', nls); 

        end
        
    end

    methods (Test)

        function testChen2020(test, useSimplifiedDiffusionPE, useSimplifiedDiffusionNE)

            doSingleTimestep = true;
            testchen2020(test, doSingleTimestep, useSimplifiedDiffusionPE, useSimplifiedDiffusionNE);
            
        end

        function testChen2020Verification(test)

            doSingleTimestep = false;
            simpleDiffusionForPE = false;
            simpleDiffusionForNE = false;
            states = testchen2020(test, doSingleTimestep, simpleDiffusionForPE, simpleDiffusionForNE);

            % Compare with PyBamm
            loadChenPybammSolution;
            [t1, u1] = deal(t, u);

            ind = cellfun(@(x) not(isempty(x)), states); 
            states = states(ind);
            time = cellfun(@(x) x.time, states);
            time = time / hour;
            Enew = cellfun(@(x) x.Control.E, states);

            x = linspace(time(1), t1(end), 100);
            pb = interp1(t1, u1, x);
            battmo = interp1(time, Enew, x);

            assert(norm(pb - battmo) < 0.077);
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
