classdef BatteryImpedanceSolver < ImpedanceSolver

    methods

        function impsolv = BatteryImpedanceSolver(inputparams, options, extrastructs)

            impsolv = impsolv@ImpedanceSolver(inputparams, options, extrastructs);
            
        end

        function setupVarNames(impsolv)
            
            impsolv.Ivarname = {'Control', 'I'};
            impsolv.Uvarname = {'Control', 'E'};
            
        end

        function jac = getJacobian(impsolv, dt)

            state = impsolv.state;
            model = impsolv.model;
            indI  = impsolv.indI;
            inds  = impsolv.inds;
            indUs = impsolv.indUs;
            
            ctrl = 'Control';

            state.(ctrl).I = 0; %steady state
            
            state0 = state; % perturbation around equilibrium, state0 is an equilibrium

            state = model.initStateAD(state);
            
            inputparams = ControlModelInputParams([]);
            model.Control = ControlModel(inputparams);
            model.Control.controlPolicy = 'None';

            drivingForces = model.getValidDrivingForces();

            % same code as BaseModel.getEquations
            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end

            state = model.applyScaling(state);
            
            for ieq = 1 : numel(model.equationVarNames)
                eqs{ieq} = model.getProp(state, model.equationVarNames{ieq});
            end
            eqs = vertcat(eqs{:});
            
            if isempty(impsolv.b)
                impsolv.b = eqs.jac{indI};
            end

            eqs.jac = eqs.jac(inds); 
            eqs = combineEquations(eqs);
            
            jac = eqs.jac{1};
            
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

