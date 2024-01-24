classdef TestChen2020 < matlab.unittest.TestCase

    methods

        function states = runChen2020(test)

        % load MRST modules
            mrstModule add ad-core mrst-gui mpfa

            % Useful abbreviations
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';

            % We create an instance of BatteryInputParams. This class is used to initiate the battery simulator and it propagates
            % all the parameters through out the submodels.

            % The input parameters can be given in json format. The json file is read and used to populate the inputparams object.
            jsonstruct_material = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));
            jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));

            jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                                           jsonstruct_geometry});

            inputparams = BatteryInputParams(jsonstruct);

            %% We setup the battery geometry ("bare" battery with no current collector).

            inputparams = setupBatteryGridFromJson(inputparams, jsonstruct);

            %%  The Battery model is initialized by sending inputparams to the Battery class constructor

            model = Battery(inputparams);

            %% We fix the input current to 5A

            model.Control.Imax = 5*ampere;

            %% We setup the schedule

            % We use different time step for the activation phase (small time steps) and the following discharging phase
            % We start with rampup time steps to go through the activation phase

            fac   = 2;
            total = 1.4*hour;
            n     = 100;
            dt0   = total*1e-6;
            times = getTimeSteps(dt0, n, total, fac);
            dt    = diff(times);
            dt    = dt(1 : end);
            step  = struct('val', dt, 'control', ones(size(dt)));

            % we setup the control for the schedule
            control = model.Control.setupScheduleControl();

            % This control is used to set up the schedule
            schedule = struct('control', control, 'step', step);

            %%  We setup the initial state

            c_ne = 29.866*mol/litre; % initial concentration at negative electrode
            c_pe = 17.038*mol/litre; % initial concentration at positive electrode
            initstate = initStateChen2020(model, c_ne, c_pe);

            %% We setup the solver
            % Default values would typically work.

            nls = NonLinearSolver();
            % Change default maximum iteration number in nonlinear solver
            nls.maxIterations = 10;
            % Change default behavior of nonlinear solver, in case of error
            nls.errorOnFailure = false;
            linearsolver = 'direct';
            switch linearsolver
              case 'agmg'
                mrstModule add agmg
                nls.LinearSolver = AGMGSolverAD('verbose', true, 'reduceToCell', true);
                nls.LinearSolver.tolerance = 1e-3;
                nls.LinearSolver.maxIterations = 30;
                nls.maxIterations = 10;
                nls.verbose = 10;
              case 'battery'
                nls.LinearSolver = LinearSolverBatteryExtra('verbose', false, 'reduceToCell', false, 'verbosity', 3, 'reuse_setup', false, 'method', 'direct');
                nls.LinearSolver.tolerance=0.5e-4*2;
              case 'direct'
                % nothing to setup in this case
              otherwise
                error('linear solver not recognized.');
            end

            % Change default tolerance for nonlinear solver
            model.nonlinearTolerance = 1e-5;
            % Set verbosity
            model.verbose = false;

            %% We run simulation

            [~, states, ~] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

        end

    end

    methods (Test)

        function testChen2020(test)

            states = test.runChen2020();

            % Load PyBamm data
            loadChenPybammSolution;

            [t1, u1] = deal(t, u);

            ind = cellfun(@(x) not(isempty(x)), states);
            states = states(ind);
            time = cellfun(@(x) x.time, states);
            time = time / hour;
            Enew = cellfun(@(x) x.Control.E, states);

            x = linspace(time(1), t1(end), 1000);
            pb = interp1(t1, u1, x);
            battmo = interp1(time, Enew, x);

            assert(norm(pb - battmo) / norm(pb) < 0.00285);

        end

    end

end


%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
