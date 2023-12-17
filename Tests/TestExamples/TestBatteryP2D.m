classdef TestBatteryP2D < matlab.unittest.TestCase

    properties (TestParameter)

        controlPolicy              = {'CCCV', 'CCDischarge'};
        use_thermal                = {true, false};
        include_current_collectors = {true, false};
        diffusionModelType         = {'simple', 'full'};
        testSize                   = {'short', 'long'};
        createReferenceData        = {false};
        compareWithReferenceData   = {true};

    end

    methods

        function states = test1d(test, controlPolicy, use_thermal, include_current_collectors, diffusionModelType, testSize, varargin)

            mrstModule add ad-core mrst-gui mpfa
            
            jsonfile = fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json');
            json = parseBattmoJson(jsonfile);

            % Change json params
            params = {'include_current_collectors', include_current_collectors, ...
                      'use_thermal', use_thermal};
            params = [params, {'NegativeElectrode.Coating.ActiveMaterial.diffusionModelType', diffusionModelType}];
            params = [params, {'PositiveElectrode.Coating.ActiveMaterial.diffusionModelType', diffusionModelType}];

            % Validation doesn't run in parallel (probably due to
            % writing to file). It doesn't run if python is not
            % available.
            try
                serial = isempty(getCurrentTask());
            catch
                serial = true;
            end
            has_python = pyenv().Version ~= "";
            validate = serial & has_python;
            validate = false;
            json = updateJson(json, params, 'validate', validate);

            inputparams = BatteryInputParams(json);

            use_cccv = strcmpi(controlPolicy, 'CCCV');
            if use_cccv
                cccvstruct = struct( 'controlPolicy'     , 'CCCV'       ,  ...
                                     'initialControl'    , 'discharging', ...
                                     'CRate'             , 1            , ...
                                     'lowerCutoffVoltage', 2.4          , ...
                                     'upperCutoffVoltage', 4.1          , ...
                                     'dIdtLimit'         , 0.01         , ...
                                     'dEdtLimit'         , 0.01);
                cccvinputparams = CcCvControlModelInputParams(cccvstruct);
                inputparams.Control = cccvinputparams;
            end

            % We define some shorthand names for simplicity.
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            elyte   = 'Electrolyte';
            thermal = 'ThermalModel';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            sd      = 'SolidDiffusion';
            ctrl    = 'Control';

            %% Setup the geometry and computational grid
            % Here, we setup the 1D computational grid that will be used for the
            % simulation. The required discretization parameters are already included
            % in the class BatteryGeneratorP2D.
            gen = BatteryGeneratorP2D();

            % Now, we update the inputparams with the properties of the grid.
            inputparams = gen.updateBatteryInputParams(inputparams);

            %%  Initialize the battery model.
            % The battery model is initialized by sending inputparams to the Battery class
            % constructor. see :class:`Battery <Battery.Battery>`.
            model = Battery(inputparams);
            model.AutoDiffBackend = AutoDiffBackend();

            %% Setup schedule
            
            step    = model.Control.setupScheduleStep();
            control = model.Control.setupScheduleControl();
            
            % This control is used to set up the schedule
            schedule = struct('control', control, 'step', step);

            %% Setup the initial state of the model
            % The initial state of the model is dispatched using the
            % model.setupInitialState()method.
            initstate = model.setupInitialState();

            %% Setup the properties of the nonlinear solver
            nls = NonLinearSolver();
            % Change default maximum iteration number in nonlinear solver
            nls.maxIterations = 10;
            % Change default behavior of nonlinear solver, in case of error
            nls.errorOnFailure = false;
            nls.timeStepSelector=StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);
            % Change default tolerance for nonlinear solver
            model.nonlinearTolerance = 1e-3*model.Control.Imax;
            % Set verbosity
            model.verbose = true;

            %% Run the simulation
            [~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

        end

    end

    methods (Test)

        function testBattery(test                      , ...
                             controlPolicy             , ...
                             use_thermal               , ...
                             include_current_collectors, ...
                             diffusionModelType        , ...
                             testSize                  , ...
                             createReferenceData       , ...
                             compareWithReferenceData)


            states = test1d(test, controlPolicy, use_thermal, include_current_collectors, diffusionModelType, testSize);

            filename = sprintf('TestBatteryP2D-%s-%d-%d-%s-%s.mat', controlPolicy, use_thermal, include_current_collectors, diffusionModelType, testSize);
            filename = fullfile(battmoDir(), 'Tests', 'TestExamples', 'ReferenceData', filename);

            if createReferenceData
                refstate = states{end};
                save(filename, 'refstate');
            elseif compareWithReferenceData
                load(filename);
                verifyStruct(test, states{end}, refstate);
            else
                % do nothing
            end

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
