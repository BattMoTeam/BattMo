classdef ElectrolyserImpedanceSolver < ImpedanceSolver

    methods

        function impsolv = ElectrolyserImpedanceSolver(inputparams, options, extrastructs)

            impsolv = impsolv@ImpedanceSolver(inputparams, options, extrastructs);
            
            % other choice for state initializaztion is 'given state'
            options = setDefaultJsonStructField(options, {'stateInitialization', 'initializationSetup'}, 'given current');

            if strcmp(getJsonStructField(options, {'stateInitialization', 'initializationSetup'}), 'given current')
                options = setDefaultJsonStructField(options, {'stateInitialization', 'computeSteadyState'}, true);
            end
            
        end

        function setupVarNames(impsolv)
            
            oer = 'OxygenEvolutionElectrode';
            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            impsolv.Ivarname = {oer, ctl, 'I'};
            impsolv.Uvarname = {oer, ptl, 'E'};
            
        end


        function drivingForces = setupDrivingForces(impsolv)
            
            drivingForces = impsolv.model.getValidDrivingForces();
                      
        end
        
        function state = setupIvalue(impsolv, state)
            
            %% TO BE UPDATED

            options = impsolv.options

            if strcmp(getStructField(options, {'stateInitialization', 'initializationSetup'}), 'given current')
                
                I = getStructField(options, {'stateInitialization', 'current'});

                assert(isAssigned(I), 'A value for the current must be given');
                
                state = model.setProp(state, impsolv.Ivarname, I);

            end

        end

        function setupSteadyState(impsolv)

            ctrl = 'Control';

            inputparams = impsolv.inputparams;
            options     = impsolv.options;
            
            inputparams.(ctrl) = CCDischargeControlModelInputParams([]);;
            inputparams.(ctrl).lowerCutoffVoltage = 3; % not used but needed for proper initialization

            initsetup = getJsonStructField(options, {'stateInitialization', 'initializationSetup'});

            switch initsetup

              case 'soc'
                
                inputparams.SOC = options.stateInitialization.soc;
                model = GenericBattery(inputparams);
                state = model.setupInitialState();

              case 'given state'

                state = impsolv.extrastructs.initstate;
                model = GenericBattery(inputparams);
                
              otherwise
                
                error('initsetup not found');
                
            end

            if getJsonStructField(options, {'stateInitialization', 'computeSteadyState'});
                
                state.(ctrl).I = 0;
                
                N         = getJsonStructField(impsolv.options, {'stateInitialization', 'numberOfTimeSteps'});
                totalTime = 10*hour;
                nr        = getJsonStructField(impsolv.options, {'stateInitialization', 'numberOfRampupSteps'});
                
                dt = rampupTimesteps(totalTime, totalTime/N, nr);

                step.val = dt;
                step.control = ones(numel(dt), 1);

                control.src = @(time, Imax) 0;

                schedule = struct('step'   , step, ...
                                  'control', control);
                
                [~, states, report] = simulateScheduleAD(state, model, schedule);

                impsolv.state = states{end};

            else

                state.time = 0; % value appears to be needed
                impsolv.state = state;
                
            end

        end
        
    end
    
end

