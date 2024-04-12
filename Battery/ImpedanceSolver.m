classdef ImpedanceSolver

    properties (SetAccess = private)

        model
        state

        %% helpers quantity (see compute Impedance to see how they are used)
        %
        DM
        DA
        b
        indUs
        
    end

    methods

        function impsolv = ImpedanceSolver(model, state, varargin)

            opt = struct('setupSteadyState', false)
            opt = merge_options(opt, varargin{:});

            impsolv.model = model;

            if nargin < 2
                state = model.setupInitialState();
                fprintf('We adjust state to equilibrim\n');
                opt.setupSteadyState = true;
            end

            if opt.setupSteadyState
                impsolv.state = impsolv.setupSteadyState(impsolv, state);
            end

            impsolv = setupHelpers();
            
        end

        function Z = computeImpedance(impsolv, omegas)
            
            DM    = impsolv.DM;
            DA    = impsolv.DA;
            b     = impsolv.b;
            indUs = impsolv.indUs;

            for iomega = 1 : numel(omegas)

                omega = omegas(iomega);

                A = (i*omega*DM + DA);
                x = A\b;

                Z = x(indUs(1) : indUs(2));
                
            end
            
        end
        
        function impsolv = setupHelpers()

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
        
        function state = setupSteadyState(impsolv, state)

            model = moldel.impsolv;
            
            ctrlinputparams = CCDischargeControlModelInputParams();
            model.Control = CCDischargeControlModel(ctrlinputparams);
            
            % This is "fragile" operation. In general it is not wise to change submodel directly, as we may miss some
            % consistency (this change is not supported).
            model.Control.Imax = 0;

            N = 10;
            totalTime = 10*hour;
            dt = rampupTimesteps(totalTime, totalTime/N, 3);

            step.val = dt;
            step.control = ones(numel(dt), 1);

            control = model.Control.setupScheduleControl();

            schedule = struct('step'   , step, ...
                              'control', control);
            
            [~, states, report] = simulateScheduleAD(state, model, schedule);

            state = states{end};
            
        end
        
    end
end
