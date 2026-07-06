classdef ImpedanceSolver < handle

    properties (SetAccess = protected)

        model
        inputparams
        state

        I
        
        %% options
        options
        extrastructs

        Ivarname
        Uvarname
        eqIvarname
        
        %% helpers quantity (see compute Impedance to see how they are used)
        %
        DM
        DA
        b
        
        indUs
        indI
        indEqI
        indEqs
        indPvs
        
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

            if strcmp(getStructField(options, {'stateInitialization', 'initializationSetup'}), 'given state')
                options = setDefaultStructField(options, {'stateInitialization', 'computeSteadyState'}, true);
            end
            
            % options used if computing steady state
            options = setDefaultStructField(options, {'stateInitialization', 'numberOfRampupSteps'}, 3); 
            options = setDefaultStructField(options, {'stateInitialization', 'numberOfTimeSteps'}  , 10);
            
            impsolv.options = options;

            impsolv.extrastructs = extrastructs;
            impsolv.inputparams  = inputparams;

            impsolv.setupVarNames();
            impsolv.setupModel();
            impsolv.setupIvalue();
            impsolv.setupSteadyState();
            impsolv.setupHelpers();
            
        end

        function setupModel(impsolv)
        %% virtual function
        end

        function setupVarNames(impsolv)
        %% virtual function
        end

        function state = setupIvalue(impsolv, state)
        %% Virtual function
        end

        function drivingForces = setupDrivingForces(impsolv)
        %% virtual function            
        end

        function setupSteadyState(impsolv)
        %% virtual function
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

            impsolv.setupIndices();
            impsolv.setupMatrices();
            
        end
        
        function setupMatrices(impsolv)
            
            dt = impsolv.options.dt;

            jac1 = impsolv.getJacobian(dt);
            jac2 = impsolv.getJacobian(2*dt);
            
            impsolv.DM = 2*dt*(jac1 - jac2);
            impsolv.DA = 2*jac2 - jac1;

        end

        function jac = getJacobian(impsolv, dt)

            % is called first because, may modify model
            drivingForces = impsolv.setupDrivingForces();
            
            state = impsolv.state;
            model = impsolv.model;
            
            state = impsolv.model.setProp(state, impsolv.Ivarname, impsolv.I);
            
            state0 = state; % perturbation around equilibrium, state0 is an equilibrium

            state = model.initStateAD(state);

            % same code as BaseModel.getEquations
            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end

            state = model.applyScaling(state);
            
            for ieq = 1 : numel(model.equationVarNames)
                eqs{ieq} = model.getProp(state, model.equationVarNames{ieq});
            end
            eqs = eqs(impsolv.indEqs);
            eqs = vertcat(eqs{:});
            
            if isempty(impsolv.b)
                impsolv.b = eqs.jac{impsolv.indI};
            end

            eqs.jac = eqs.jac(impsolv.indPvs); 
            eqs = combineEquations(eqs);
            
            jac = eqs.jac{1};
            
        end

        
        function setupIndices(impsolv)

            model = impsolv.model;
            
            Ivarname   = impsolv.Ivarname;
            Uvarname   = impsolv.Uvarname;
            eqIvarname = impsolv.eqIvarname;
            
            indI = model.getIndexPrimaryVariable(Ivarname);
            impsolv.indI = indI;
            
            indEqI = model.getIndexEquationVarName(eqIvarname);
            impsolv.indEqI = indEqI;

            indPvs = true(numel(model.getPrimaryVariableNames), 1);
            indPvs(indI) = false;
            indPvs = find(indPvs);

            impsolv.indPvs = indPvs;

            indEqs = true(numel(model.getPrimaryVariableNames), 1);
            indEqs(indEqI) = false;
            indEqs = find(indEqs);

            impsolv.indEqs = indEqs;
            
            state = model.initStateAD(impsolv.state); % just needed to get AD sample
            adsample = model.getProp(state, Uvarname);
            indUs = model.getRangePrimaryVariable(adsample, Uvarname);

            impsolv.indUs = indUs;
            
        end
        
    end
    
end

