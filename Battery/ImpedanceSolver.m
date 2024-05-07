classdef ImpedanceSolver

    properties (SetAccess = private)

        model
        inputparams
        state

        %% options
        options
        
        %% helpers quantity (see compute Impedance to see how they are used)
        %
        DM
        DA
        b
        indUs
        
    end

    methods

        function impsolv = ImpedanceSolver(inputparams, varargin)

            default_options = struct('computeSteadyState' , true, ...
                                     'soc'                , []  , ...
                                     'initstate'          , []  , ...
                                     'numberOfRampupSteps', 3   , ...
                                     'numberOfTimeSteps'  , 10);

            impsolv.options = merge_options(default_options, varargin{:});

            ctrl = 'Control';

            impsolv.inputparams = inputparams;

            inputparams.(ctrl) = ImpedanceControlModelInputParams([]);
            impsolv.model = ImpedanceBattery(inputparams);

            impsolv = impsolv.setupSteadyState();
            impsolv = impsolv.setupHelpers();
            
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
        
        function impsolv = setupHelpers(impsolv)

            state = impsolv.state;
            model = impsolv.model;
            
            ctrl = 'Control';

            state.(ctrl).omega = 1;
            state.(ctrl).I = 0;
            
            state = model.initStateAD(state);

            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end

            state = model.applyScaling(state);
            
            for ieq = 1 : numel(model.equationVarNames)
                eqs{ieq} = model.getProp(state, model.equationVarNames{ieq});
            end
            eqs = vertcat(eqs{:});
            
            indI = model.getIndexPrimaryVariable({ctrl, 'I'});

            inds = true(numel(model.getPrimaryVariableNames), 1);
            inds(indI) = false;
            inds = find(inds);

            b = eqs.jac{indI};

            eqs.jac = eqs.jac(inds); 
            eqs = combineEquations(eqs);

            indUs = model.getRangePrimaryVariable(state.(ctrl).E, {ctrl, 'E'});
            
            A = eqs.jac{1};

            impsolv.DM    = imag(A);
            impsolv.DA    = real(A);
            impsolv.b     = b;
            impsolv.indUs = indUs;
            
        end

        function impsolv = setupSteadyState(impsolv)

            ctrl = 'Control';

            inputparams = impsolv.inputparams;
            options     = impsolv.options;
            
            inputparams.(ctrl) = CCDischargeControlModelInputParams([]);;
            inputparams.(ctrl).lowerCutoffVoltage = 3; % not used but needed for proper initialization

            initcase = [];
            if ~isempty(options.soc)
                initcase = 'soc';
            elseif ~isempty(options.initstate)
                initcase = 'state';
            else
                error('initcase not found');
            end


            switch initcase

              case 'soc'
                
                inputparams.SOC = options.soc;
                model = Battery(inputparams);
                state = model.setupInitialState();

              case 'state'

                state = options.initstate;
                model = Battery(inputparams);
                
              otherwise
                
                error('initcase not found');
                
            end


            if options.computeSteadyState
                
                state.(ctrl).I = 0;
                
                N         = impsolv.options.numberOfTimeSteps;
                totalTime = 10*hour;
                nr        = impsolv.options.numberOfRampupSteps;
                
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
