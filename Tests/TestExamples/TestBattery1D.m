classdef TestBattery1D < matlab.unittest.TestCase

    properties (TestParameter)
        
        jsonfile                  = {fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json')};
        controlPolicy             = {'CCCV', 'IEswitch'};
        use_thermal               = {true, false};
        include_current_collector = {true, false};
        diffusionmodel            = {'none', 'simplified', 'full'};
        testSize                  = {'short', 'long'};
        createReferenceData       = {false};
    end
    
    methods

        function states = test1d(test, jsonfile, controlPolicy, use_thermal, include_current_collector, diffusionModelType, testSize, varargin)
            
            jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));


            jsonstruct.include_current_collector = include_current_collector;
            jsonstruct.use_thermal = use_thermal;

            if strcmp(diffusionModelType, 'none')
                json.use_particle_diffusion = false;
            else
                json.use_particle_diffusion = true;
                json.NegativeElectrode.ActiveMaterial.diffusionModelType = diffusionModelType;
                json.PositiveElectrode.ActiveMaterial.diffusionModelType = diffusionModelType;
            end
            
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


            % We define some shorthand names for simplicity.
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            elyte   = 'Electrolyte';
            thermal = 'ThermalModel';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            sd      = 'SolidDiffusion';
            ctrl    = 'Control';

            %% Setup the geometry and computational mesh
            % Here, we setup the 1D computational mesh that will be used for the
            % simulation. The required discretization parameters are already included
            % in the class BatteryGenerator1D. 
            gen = BatteryGenerator1D();

            % Now, we update the paramobj with the properties of the mesh. 
            paramobj = gen.updateBatteryInputParams(paramobj);

            %%  Initialize the battery model. 
            % The battery model is initialized by sending paramobj to the Battery class
            % constructor. see :class:`Battery <Battery.Battery>`.
            model = Battery(paramobj);
            model.AutoDiffBackend= AutoDiffBackend();

            %% Compute the nominal cell capacity and choose a C-Rate
            % The nominal capacity of the cell is calculated from the active materials.
            % This value is then combined with the user-defined C-Rate to set the cell
            % operational current. 

            CRate = model.Control.CRate;

            %% Setup the time step schedule 
            % Smaller time steps are used to ramp up the current from zero to its
            % operational value. Larger time steps are then used for the normal
            % operation.
            switch model.(ctrl).controlPolicy
              case 'CCCV'
                total = 3.5*hour/CRate;
              case 'IEswitch'
                total = 1.4*hour/CRate;
              otherwise
                error('control policy not recognized');
            end

            n     = 100;
            dt    = total/n;

            switch testSize
              case 'long'
                % do nothing
              case 'short'
                n = 10;
              otherwise
                error('testSize not recognized')
            end
            
            
            step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

            % we setup the control by assigning a source and stop function.
            % control = struct('CCCV', true); 
            %  !!! Change this to an entry in the JSON with better variable names !!!

            switch model.Control.controlPolicy
              case 'IEswitch'
                tup = 0.1; % rampup value for the current function, see rampupSwitchControl
                srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                            model.Control.Imax, ...
                                                            model.Control.lowerCutoffVoltage);
                % we setup the control by assigning a source and stop function.
                control = struct('src', srcfunc, 'IEswitch', true);
              case 'CCCV'
                control = struct('CCCV', true);
              otherwise
                error('control policy not recognized');
            end

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
            NLS.errorOnFailure = false; 
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

        function testBattery(test, jsonfile, controlPolicy, use_thermal, include_current_collector, diffusionmodel, testSize, createReferenceData)

            states = test1d(test, jsonfile, controlPolicy, use_thermal, include_current_collector, diffusionmodel, testSize);

            filename = sprintf('TestBattery1D-%s-%d-%d-%s-%s.mat', controlPolicy, use_thermal, include_current_collector, diffusionmodel, testSize);
            filename = fullfile(battmoDir(), 'Tests', 'TestExamples', 'ReferenceData', filename);

            if createReferenceData
                refstate = states{end};
                save(filename, 'refstate');
            else
                load(filename);
                verifyStruct(test, states{end}, refstate);
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
