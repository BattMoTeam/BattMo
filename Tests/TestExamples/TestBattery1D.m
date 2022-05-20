classdef TestBattery1D < matlab.unittest.TestCase

    properties (TestParameter)
        
        jsonfile = {'ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json'};
        controlPolicy = {'CCCV', 'IEswitch'};
        
    end
    
    methods

        function test = TestBattery1D()
            mrstModule reset
            mrstModule add ad-core mrst-gui mpfa
        end
        
        function states = test1d(test, jsonfile, controlPolicy, varargin)

            % Setup properties
            jsonstruct = parseBattmoJson(jsonfile);
            paramobj = BatteryInputParams(jsonstruct);

            use_cccv = strcmpi(controlPolicy, 'CCCV');
            if use_cccv
                cccvstruct = struct( 'controlPolicy'     , 'CCCV',  ...
                                     'CRate'             , 1         , ...
                                     'lowerCutoffVoltage', 2         , ...
                                     'upperCutoffVoltage', 4.1       , ...
                                     'dIdtLimit'         , 0.01      , ...
                                     'dEdtLimit'         , 0.01);
                cccvparamobj = CcCvControlModelInputParams(cccvstruct);
                paramobj.Control = cccvparamobj;
            end

            % Setup geometry and mesh
            gen = BatteryGenerator1D();
            paramobj = gen.updateBatteryInputParams(paramobj);

            % Initialize the model
            model = Battery(paramobj);
            model.AutoDiffBackend = AutoDiffBackend();

            % Choose C-rate
            CRate = model.Control.CRate;

            % Setup time step schedule
            switch model.Control.controlPolicy
              case 'CCCV'
                total = 3.5*hour/CRate;
              case 'IEswitch'
                total = 1.4*hour/CRate;
              otherwise
                error('control policy not recognized');
            end

            n    = 100;
            dt   = total/n;
            step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));
            
            switch model.Control.controlPolicy
              case 'IEswitch'
                tup = 0.1;
                srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                            model.Control.Imax, ...
                                                            model.Control.lowerCutoffVoltage);
                control = struct('src', srcfunc, 'IEswitch', true);
              case 'CCCV'
                control = struct('CCCV', true);
              otherwise
                error('control policy not recognized');
            end

            schedule = struct('control', control, 'step', step); 

            % Setup initial state
            initstate = model.setupInitialState();

            % Setup solver properties
            nls = NonLinearSolver(); 
            nls.errorOnFailure = true;
            nls.timeStepSelector = StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);
            model.nonlinearTolerance = 1e-3*model.Control.Imax;
            model.verbose = false;

            [~, states] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', false, 'NonLinearSolver', nls); 

        end
        
    end

    methods (Test)

        function testBattery(test, jsonfile, controlPolicy)

            states = test1d(test, jsonfile, controlPolicy);
            verifyStruct(test, states{end}, sprintf('TestBattery1D%s', controlPolicy));
            
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
