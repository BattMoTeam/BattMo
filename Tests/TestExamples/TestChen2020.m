classdef TestChen2020 < matlab.unittest.TestCase

    methods

        
        function states = testchen2020(test)
            
            jsonstruct = parseBattmoJson(fullfile('ParameterData','ParameterSets','Chen2020','chen2020_lithium_ion_battery.json'));

            paramobj = BatteryInputParams(jsonstruct);
            paramobj.include_current_collectors = false;
            paramobj = paramobj.validateInputParams();
            
            % Some shorthands used for the sub-models
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            itf   = 'Interface';
            sd    = 'SolidDiffusion';
            elyte = 'Electrolyte';
            ctrl  = 'Control';
            
            %% We setup the battery geometry ("bare" battery with no current collector).
            gen = BareBatteryGenerator3D();
            % We update pamobj with grid data
            paramobj = gen.updateBatteryInputParams(paramobj);

            %%  The Battery model is initialized by sending paramobj to the Battery class constructor 

            model = Battery(paramobj);

            %% We fix the input current to 5A

            model.Control.Imax = 5;

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

            % We set up a stopping function. Here, the simulation will stop if the output voltage reach a value smaller than 2. This
            % stopping function will not be triggered in this case as we switch to voltage control when E=3.6 (see value of inputE
            % below).

            tup = 0.1; % rampup value for the current function, see rampupSwitchControl
            srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                        model.Control.Imax, ...
                                                        model.Control.lowerCutoffVoltage);
            % we setup the control by assigning a source and stop function.
            control = struct('src', srcfunc, 'IEswitch', true);

            % This control is used to set up the schedule
            schedule = struct('control', control, 'step', step); 

            %%  We setup the initial state

            nc = model.G.cells.num;
            T = model.initT;
            initstate.ThermalModel.T = T*ones(nc, 1);

            bat = model;

            initstate = model.updateTemperature(initstate);

            % we setup negative electrode initial state
            nitf = bat.(ne).(am).(itf); 

            % We bypass the solid diffusion equation to set directly the particle surface concentration
            c = 29866.0;
            if strcmp(model.(ne).(am).diffusionModelType, 'simple')
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

            % we setup positive electrode initial state

            pitf = bat.(pe).(am).(itf); 

            c = 17038.0;

            if strcmp(model.(pe).(am).diffusionModelType, 'simple')
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

            % Setup nonlinear solver 
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
                %nls.LinearSolver = LinearSolverBatteryExtra('verbose', false, 'reduceToCell', true,'verbosity',3,'reuse_setup',false,'method','matlab_p_gs');
                nls.LinearSolver = LinearSolverBatteryExtra('verbose', false, 'reduceToCell', false,'verbosity',3,'reuse_setup',false,'method','direct');
                nls.LinearSolver.tolerance=0.5e-4*2;          
              case 'direct'
                disp('standard direct solver')
              otherwise
                error()
            end

            % Change default tolerance for nonlinear solver
            model.nonlinearTolerance = 1e-5; 
            % Set verbosity
            model.verbose = false;

            model.AutoDiffBackend= AutoDiffBackend();

            % Run simulation
            [wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

        end
        
    end

    methods (Test)

        function testChen2020(test)

            states = testchen2020(test);

            % Compare with PyBamm
            loadChenPybammSolution;
            [t1, u1] = deal(t, u);

            ind = cellfun(@(x) not(isempty(x)), states); 
            states = states(ind);
            time = cellfun(@(x) x.time, states);
            time = time / hour;
            Enew = cellfun(@(x) x.Control.E, states);

            x = linspace(time(1), t1(end), 100000);
            pb = interp1(t1, u1, x);
            battmo = interp1(time, Enew, x);

            assert(norm(pb - battmo) / norm(pb) < 0.00285);
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
