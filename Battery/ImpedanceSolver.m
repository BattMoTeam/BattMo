classdef ImpedanceSolver < handle

    properties (SetAccess = private)

        model
        inputparams
        state

        %% options
        options
        extrastructs
        
        %% helpers quantity (see compute Impedance to see how they are used)
        %
        DM
        DA
        b
        
        indUs
        indI
        inds
        
    end

    methods

        function impsolv = ImpedanceSolver(inputparams, options, extrastructs)
            
            if nargin < 3
                extrastructs = [];
                if nargin < 2
                    options = [];
                end
            end

            % reference dt used to compute jacobians
            options = setDefaultStructField(options, {'dt'}, 1e3);
            
            % other choice for state initializaztion is 'given state'
            options = setDefaultStructField(options, {'stateInitialization', 'initializationSetup'}, 'soc');  
                                                                                                                             
            % default soc value
            options = setDefaultStructField(options, {'stateInitialization', 'soc'}, 1);

            switch getStructField(options, {'stateInitialization', 'initializationSetup'})
              case 'soc'
                options = setDefaultStructField(options, {'stateInitialization', 'computeSteadyState'}, false);
              case 'given state'
                options = setDefaultStructField(options, {'stateInitialization', 'computeSteadyState'}, true);
              otherwise
                error('initializationSetup not recognized');
            end
            
            % options used if computing steady state
            options = setDefaultStructField(options, {'stateInitialization', 'numberOfRampupSteps'}, 3); 
            options = setDefaultStructField(options, {'stateInitialization', 'numberOfTimeSteps'}  , 10);
            
            impsolv.options = options;

            impsolv.extrastructs = extrastructs;

            ctrl = 'Control';

            impsolv.inputparams = inputparams;

            inputparams.(ctrl) = CCDischargeControlModelInputParams([]);
            impsolv.model = GenericBattery(inputparams);

            impsolv.setupSteadyState();
            impsolv.setupHelpers();
            
        end

        function Z = computeImpedance(impsolv, omegas)
            
            DM    = impsolv.DM;
            DA    = impsolv.DA;
            b     = impsolv.b;
            indUs = impsolv.indUs;

            for iomega = 1 : numel(omegas)

                omega = omegas(iomega);

                A = (i*2*pi*omega*DM + DA);
                x = A\b;

                Z(iomega) = x(indUs(1) : indUs(2));
                
            end
            
        end
        
        function setupHelpers(impsolv)

            model = impsolv.model;
            
            ctrl = 'Control';
            
            indI = model.getIndexPrimaryVariable({ctrl, 'I'});

            impsolv.indI = indI;
            
            inds = true(numel(model.getPrimaryVariableNames), 1);
            inds(indI) = false;
            inds = find(inds);

            impsolv.inds = inds;
            
            state = model.initStateAD(impsolv.state); % just needed to get AD sample
            indUs = model.getRangePrimaryVariable(state.(ctrl).E, {ctrl, 'E'});

            impsolv.indUs = indUs;

            dt = impsolv.options.dt;

            jac1 = impsolv.getJacobian(dt);
            jac2 = impsolv.getJacobian(2*dt);
            
            impsolv.DM = 2*dt*(jac1 - jac2);
            impsolv.DA = 2*jac2 - jac1;

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

                state.time = 0; % value appears to be needed
                impsolv.state = state;
                
            end

        end
        
    end
    
end

