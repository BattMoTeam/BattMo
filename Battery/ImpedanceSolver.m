classdef ImpedanceSolver

    properties (SetAccess = private)

        model
        inputparams
        state
        soc

        %% helpers quantity (see compute Impedance to see how they are used)
        %
        DM
        DA
        b
        indUs
        
    end

    methods

        function impsolv = ImpedanceSolver(inputparams, soc)

            ctrl = 'Control';

            impsolv.soc         = soc;
            impsolv.inputparams = inputparams;

            inputparams.(ctrl) = ImpedanceControlModelInputParams([]);
            impsolv.model = ImpedanceBattery(inputparams);

            impsolv = impsolv.setupSteadyState(soc);
            impsolv = impsolv.setupHelpers();
            
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
        
        function state = setupSteadyState(impsolv, state)

            ctrl = 'Control';

            inputparams = impsolv.inputparams;
            inputparams.(ctrl) = CCDischargeControlModelInputParams([]);;
            inputparams.(ctrl).lowerCutoffVoltage = 3; % not used but needed for proper initialization
            model = Battery(inputparams);

            if isempty(state)
                state = model.setupInitialState();
            else
                state.(ctrl).ctrlType = 'constantCurrent';
                state.(ctrl).I = 0;
            end
            
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
