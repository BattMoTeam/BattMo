classdef BatteryImpedanceSolver < ImpedanceSolver

    methods

        function impsolv = BatteryImpedanceSolver(inputparams, options, extrastructs)

            if nargin < 3
                extrastructs = [];
                if nargin < 2
                    options = [];
                end
            end

            impsolv = impsolv@ImpedanceSolver(inputparams, options, extrastructs);

            options = impsolv.options;

            % other choice for state initializaztion is 'given state'
            options = setDefaultStructField(options, {'stateInitialization', 'initializationSetup'}, 'soc');  
                                                                                                                             
            % default soc value
            options = setDefaultStructField(options, {'stateInitialization', 'soc'}, 1);
            
            if strcmp(getStructField(options, {'stateInitialization', 'initializationSetup'}), 'soc')
                options = setDefaultStructField(options, {'stateInitialization', 'computeSteadyState'}, false);
            end

        end

        function setupIvalue(impsolv)

            impsolv.I = 0;
            
        end
        
        function setupModel(impsolv)

            % We assign some control. Not used in computation
            ctrl = 'Control';
            jsonstruct.DRate = 1;
            impsolv.inputparams.(ctrl) = CCDischargeControlModelInputParams(jsonstruct);
            
            impsolv.model = GenericBattery(impsolv.inputparams);
            
        end
        
        function setupVarNames(impsolv)
            
            impsolv.Ivarname   = {'Control', 'I'};
            impsolv.Uvarname   = {'Control', 'E'};
            impsolv.eqIvarname = {'Control', 'controlEquation'};
            
        end


        function drivingForces = setupDrivingForces(impsolv)


            drivingForces = impsolv.model.getValidDrivingForces();
            schedule = impsolv.model.Control.setupSchedule();
            drivingForces.src = schedule.control.src;

        end
        
        function setupSteadyState(impsolv)

            ctrl = 'Control';

            inputparams = impsolv.inputparams;
            options     = impsolv.options;
            
            inputparams.(ctrl) = CCDischargeControlModelInputParams([]);;
            inputparams.(ctrl).lowerCutoffVoltage = 3; % not used but needed for proper initialization

            initsetup = getStructField(options, {'stateInitialization', 'initializationSetup'});

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

            if getStructField(options, {'stateInitialization', 'computeSteadyState'});
                
                state.(ctrl).I = 0;
                
                N         = getStructField(impsolv.options, {'stateInitialization', 'numberOfTimeSteps'});
                totalTime = 10*hour;
                nr        = getStructField(impsolv.options, {'stateInitialization', 'numberOfRampupSteps'});
                
                dt = rampupTimesteps(totalTime, totalTime/N, nr);

                step.val = dt;
                step.control = ones(numel(dt), 1);

                control.src = @(time, Imax) 0;

                schedule = struct('step'   , step, ...
                                  'control', control);
                
                [~, states, report] = simulateScheduleAD(state, model, schedule);

                impsolv.state = states{end};

            else

                impsolv.state = state;
                
            end

            impsolv.state.time = 0;
            
        end
        
    end
    
end

