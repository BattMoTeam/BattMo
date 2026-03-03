classdef EquilibriumConcentrationSolver < BaseModel

    properties

        % Each electrode has the following fields
        % - numberOfActiveMaterial   : number of active materials (nam)
        % - volumes                  : vector (dimension nam) with volume of each active material
        % - saturationConcentrations : vector (dimension nam) with saturation concentrations for each active material
        % - thetaMaxs                : vector (dimension nam) with maximum concentation for each active material
        % - thetaMins                : vector (dimension nam) with minimum concentation for each active material
        % - computeOCPs              : cell array where each cell is a function computing the OCP from a given concentration.
        
        NegativeElectrode
        PositiveElectrode

        voltage
        totalAmount % in mol

        % Temperature (may be needed in OCP computations)
        T
        
    end


    methods

        function model = EquilibriumConcentrationSolver(batterymodel)

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';

            model.subModelNameList = {ne, pe};

            model = model.setupFromBatteryModel(batterymodel);
            model = model.equipModelForComputation();
            
        end

        function model = registerVarAndPropfuncNames(model)

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                nam = model.(elde).numberOfActiveMaterial;
                
                varnames = {};

                varnames{end + 1} = VarName({elde}, 'stoichiometries', nam);
                varnames{end + 1} = VarName({elde}, 'concentrations', nam);
                varnames{end + 1} = VarName({elde}, 'ocps', nam);
                if nam > 1
                    varnames{end + 1} = VarName({elde}, 'potentialEquations', nam - 1);
                end

                for ivarname = 1 : numel(varnames)
                    varnames{ivarname}.useCell = true;
                end
                
                model = model.registerVarNames(varnames);

                model = model.registerVarName({elde, 'amount'});
                model = model.setAsExtraVarName({elde, 'amount'});
                
            end

            varnames = {};

            varnames{end + 1} = 'massConsEq';
            varnames{end + 1} = 'totalAmount';
            varnames{end + 1} = 'voltage';
            varnames{end + 1} = 'voltageEquation';

            model = model.registerVarNames(varnames);
            
            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                nam = model.(elde).numberOfActiveMaterial;

                fn = @EquilibriumConcentrationSolver.updateConcentrations;
                inputvarnames  = {VarName({elde}, 'stoichiometries', nam)};
                outputvarnames = VarName({elde}, 'concentrations', nam);
                model = model.registerPropFunction({outputvarnames, fn, inputvarnames});

                fn = @EquilibriumConcentrationSolver.updateOCPs;
                inputvarnames  = {VarName({elde}, 'concentrations', nam)};
                outputvarnames = VarName({elde}, 'ocps', nam);
                model = model.registerPropFunction({outputvarnames, fn, inputvarnames});

                if nam > 1
                    fn = @EquilibriumConcentrationSolver.updatePotentialEquation;
                    inputvarnames  = {VarName({elde}, 'ocps', nam)};
                    outputvarnames = VarName({elde}, 'potentialEquations', nam - 1);
                    model = model.registerPropFunction({outputvarnames, fn, inputvarnames});
                end

                fn = @EquilibriumConcentrationSolver.updateAmount;
                inputvarnames  = {VarName({elde}, 'concentrations', nam)};
                outputvarnames = {elde, 'amount'};
                model = model.registerPropFunction({outputvarnames, fn, inputvarnames});
                
                
            end

            nam_ne = model.(ne).numberOfActiveMaterial;
            nam_pe = model.(pe).numberOfActiveMaterial;
            
            fn = @ConcentrationSolve.updateVoltage;
            inputnames = {VarName({ne}, 'ocps', nam_ne), ...
                          VarName({pe}, 'ocps', nam_pe)};
            model = model.registerPropFunction({'voltage', fn, inputnames});
            
            fn = @ConcentrationSolve.updateVoltageEquation;
            inputnames = {'voltage'};
            model = model.registerPropFunction({'voltageEquation', fn, inputnames});

            fn = @ConcentrationSolve.updateTotalAmount;
            inputnames = {VarName({ne}, 'concentrations', nam_ne), ...
                          VarName({pe}, 'concentrations', nam_pe)};
            model = model.registerPropFunction({'totalAmount', fn, inputnames});

            fn = @ConcentrationSolve.updateMassConsEq;
            inputnames = {'totalAmount'};
            model = model.registerPropFunction({'massConsEq', fn, inputnames});

        end

        function model = setupFromBatteryModel(model, batterymodel)
            
            model.T = batterymodel.initT;
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            itf = 'Interface';
            
            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                switch elde
                  case ne
                    gsmax = 'guestStoichiometry100';
                    gsmin = 'guestStoichiometry0';
                  case pe
                    gsmax = 'guestStoichiometry0';
                    gsmin = 'guestStoichiometry100';
                  otherwise
                    error('Electrode not recognised');
                end
                
                % For the moment we support only 2 active materials
                if ~isempty(batterymodel.(elde).(co).ActiveMaterial1)
                    
                    comodel = batterymodel.(elde).(co);
                    nam = 2;

                    satConcs    = nan(nam, 1);
                    thetaMaxs   = nan(nam, 1);
                    thetaMins   = nan(nam, 1);
                    vfs         = nan(nam, 1);
                    vs          = nan(nam, 1);
                    computeOCPs = cell(nam, 1);
                    
                    for iam = 1 : nam

                        switch iam
                          case 1
                            am = 'ActiveMaterial1';
                          case 2
                            am = 'ActiveMaterial2';
                          otherwise
                            error('iam index not accepted');
                        end
                        
                        indam = comodel.compInds.(am);

                        satConcs(iam)    = comodel.(am).(itf).saturationConcentration;
                        thetaMaxs(iam)   = comodel.(am).(itf).(gsmax);
                        thetaMins(iam)   = comodel.(am).(itf).(gsmin);
                        vfs(iam)         = comodel.volumeFractions(indam)*comodel.volumeFraction;
                        vs(iam)          = sum(comodel.G.getVolumes()*vfs(iam));
                        computeOCPs{iam} = @(c) comodel.(am).(itf).computeOCP(c/satConcs(iam));

                    end

                    model.(elde).numberOfActiveMaterial   = nam;
                    model.(elde).saturationConcentrations = satConcs;
                    model.(elde).thetaMaxs                = thetaMaxs;
                    model.(elde).thetaMins                = thetaMins;
                    model.(elde).volumes                  = vs;
                    model.(elde).computeOCPs              = computeOCPs;
                    
                else

                    comodel = batterymodel.(elde).(co);
                    
                    am = 'ActiveMaterial';

                    indam = comodel.compInds.(am);

                    satConcs  = comodel.(am).(itf).saturationConcentration;
                    thetaMaxs = comodel.(am).(itf).(gsmax);
                    thetaMins = comodel.(am).(itf).(gsmin);
                    vfs       = comodel.volumeFractions(indam)*comodel.volumeFraction;
                    vs        = sum(comodel.G.getVolumes()*vfs);

                    computeOCPs = cell(1, 1);
                    computeOCPs{1} = @(c) comodel.(am).(itf).computeOCP(c./satConcs);

                    model.(elde).numberOfActiveMaterial   = 1;
                    model.(elde).saturationConcentrations = satConcs;
                    model.(elde).thetaMaxs                = thetaMaxs;
                    model.(elde).thetaMins                = thetaMins;
                    model.(elde).volumes                  = vs;
                    model.(elde).computeOCPs              = computeOCPs;
                    
                end
                
            end
            
        end
        
        function [initstate, model] = setupInitialState(model)
            
        % We use a fully discharge charge battery

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';

            if ~isempty(model.voltage)
                nc = numel(model.voltage);
            else
                nc = 1;
            end
            
            nam = model.(ne).numberOfActiveMaterial;
            for iam = 1 : nam
                inistate.(ne).stoichiometries{iam} = model.(ne).thetaMins(iam)*ones(nc, 1);
            end

            nam = model.(pe).numberOfActiveMaterial;
            for iam = 1 : nam
                inistate.(pe).stoichiometries{iam} = model.(pe).thetaMaxs(iam)*ones(nc, 1);
            end

            initstate = model.evalVarName(inistate, {'totalAmount'});
            
            model.totalAmount = initstate.totalAmount;
            
        end
        
        
        

        function state = updateConcentrations(model, state)

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                nam = model.(elde).numberOfActiveMaterial;

                for iam = 1 : nam
                    state.(elde).concentrations{iam} = model.(elde).saturationConcentrations(iam).*state.(elde).stoichiometries{iam};
                end
                
            end

        end

        function state = updateAmount(model, state)

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                amount = 0*state.(elde).concentrations{1}; % dummy AD initialization
                
                elde = eldes{ielde};

                nam = model.(elde).numberOfActiveMaterial;

                concs = state.(elde).concentrations;
                
                nam = model.(elde).numberOfActiveMaterial;

                for iam = 1 : nam
                    amount = concs{iam}.*model.(elde).volumes(iam);
                end

                state.(elde).amount = amount;

            end        

        end            
        
        function state = updateOCPs(model, state)

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                nam = model.(elde).numberOfActiveMaterial;

                concs = state.(elde).concentrations;
                
                nam = model.(elde).numberOfActiveMaterial;

                for iam = 1 : nam
                    
                    ocpfunc = model.(elde).computeOCPs{iam};
                    ocps{iam} = ocpfunc(concs{iam});
                    
                end

                state.(elde).ocps = ocps;

            end        

        end

        
        function state = updatePotentialEquation(model, state)

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                nam = model.(elde).numberOfActiveMaterial;

                ocps = state.(elde).ocps;

                for iam = 1 : (nam - 1)
                    
                    state.(elde).potentialEquations{iam} = ocps{iam + 1} - ocps{1};
                    
                end

            end        
            
        end

        function state = updateVoltage(model, state)
            
            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';

            state.voltage = state.(pe).ocps{1} - state.(ne).ocps{1};
            
        end

        function state = updateVoltageEquation(model, state)

            state.voltageEquation = state.voltage - model.voltage;

        end
        
        function state = updateTotalAmount(model, state)
            
            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';

            totalAmount = 0*state.(ne).concentrations{1}; % dummy AD initialization

            eldes = {ne, pe};
            
            for ielde = 1 : numel(eldes)
                
                elde = eldes{ielde};

                nam = model.(elde).numberOfActiveMaterial;

                for iam = 1 : nam
                    totalAmount = totalAmount + state.(elde).concentrations{iam}*model.(elde).volumes(iam);
                end
                
            end
            
            state.totalAmount = totalAmount;

        end
        
        function state = updateMassConsEq(model, state)
            
            state.massConsEq = state.totalAmount - model.totalAmount;

        end

        function cleanState = addStaticVariables(model, cleanState, state)
        % nothing to do here
        end

        function [state, report] = updateState(model, state, problem, dx, forces)

            [state, report] = updateState@BaseModel(model, state, problem, dx, forces);
            
            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';

            eldes = {ne, pe};
            
            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};
                for iam = 1 : model.(elde).numberOfActiveMaterial
                    state = model.capProperty(state, {elde, 'stoichiometries', iam}, model.(elde).thetaMins(iam), model.(elde).thetaMaxs(iam));
                end

            end
            
        end

        function plotOCPs(model, fig)

            if nargin > 1
                figure(fig)
            else
                figure
            end
            
            hold on

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';

            eldes = {ne, pe};

            N = 100;

            theta = linspace(0, 1, N);

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};
                
                nam = model.(elde).numberOfActiveMaterial;

                theta_elde = {};
                
                for iam = 1 : nam

                    theta_elde{iam} = linspace(model.(elde).thetaMins(iam), model.(elde).thetaMaxs(iam), N);
                    
                    state1.(elde).stoichiometries{iam} = theta;
                    state2.(elde).stoichiometries{iam} = theta_elde{iam};
                    
                end

                thetas.(elde) = theta_elde;
                
            end

            state1 = model.evalVarName(state1, 'ocps');
            state2 = model.evalVarName(state2, 'ocps');
            
            
            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};
                nam = model.(elde).numberOfActiveMaterial;
                
                for iam = 1 : nam

                    plot(theta, state1.(elde).ocps{iam});
                    plot(thetas.(elde){iam}, state2.(elde).ocps{iam}, 'linewidth', 3);
                    
                end

            end
            
            
        end
        
        function [state, failure, model] = computeConcentrations(model, voltage, varargin)

            opt = struct('verbose', false);
            opt = merge_options(opt, varargin{:});


            model.voltage = voltage;
            [initstate, model] = setupInitialState(model);

            nls = NonLinearSolver();

            if opt.verbose
                model.verbose = true;
                nls.continueOnFailure = false;
                nls.verbose = true;
            end
            
            [state, failure, report] = nls.solveMinistep(model, initstate, initstate, [], []);

        end

        function [funcs, options, extras] = setupIpOpt(model, optimcase)

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';

            eldes = {ne, pe};

            function x = convertFromState(state)
                nam1 = 1;
                for ielde = 1 : numel(eldes)
                    elde = eldes{ielde};
                    nam2 = nam1 + model.(elde).numberOfActiveMaterial - 1;
                    for iam = nam1 : nam2
                        x(iam) = state.(elde).stoichiometries{iam - nam1 + 1};
                    end
                    nam1 = nam2 + 1;
                end
            end

            function state = convertToState(x)
                nam1 = 1;
                for ielde = 1 : numel(eldes)
                    elde = eldes{ielde};
                    nam2 = nam1 + model.(elde).numberOfActiveMaterial - 1;
                    for iam = nam1 : nam2
                        state.(elde).stoichiometries{iam - nam1 + 1}= x(iam);
                    end
                    nam1 = nam2 + 1;
                end
            end

            extras.convertToState   = @(x) convertToState(x);
            extras.convertFromState = @(state) convertFromState(state);

            switch optimcase
              case 'max'
                coef = -1;
              case 'min'
                coef = 1;
              otherwise
                error('optimcase not recognized');
            end
            
            function v = computeVoltage(x)

                state = convertToState(x);
                state = model.evalVarName(state, 'voltage');
                % We change sign as ipopt computes the minimum
                v = coef*state.voltage;
                
            end


            function derv = computeDerVoltage(x)

                x = initVariablesADI(x);
                v = computeVoltage(x);
                derv = full(v.jac{1});
                
            end

            % Define the call list to update the constraint equations
            
            cg = model.computationalGraph;
            
            funcCallList = {};
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                if model.(elde).numberOfActiveMaterial > 1
                    funcCallList = horzcat(funcCallList, ...
                                           cg.getPropFunctionCallList({elde, 'potentialEquations'}));
                end
            end
            
            funcCallList = horzcat(funcCallList, ...
                                   cg.getPropFunctionCallList({'massConsEq'}));

            funcCallList = unique(funcCallList, 'stable');

            function v = computeConstraints(x)

                state = convertToState(x);
                for ifunc = 1 : numel(funcCallList)
                    eval(funcCallList{ifunc});
                end

                v = state.massConsEq;
                
                for ielde = 1 : numel(eldes)
                    elde = eldes{ielde};
                    if model.(elde).numberOfActiveMaterial > 1
                        v = vertcat(v, ...
                                    vertcat(state.(elde).potentialEquations{:}));
                    end
                end
                
            end

            function derv = computeDerConstraints(x)

                x = initVariablesADI(x);
                v = computeConstraints(x);
                derv = v.jac{1};

            end

            nam = 0;
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                nam = nam + model.(elde).numberOfActiveMaterial;
            end

            function v = getConstraintJacobianStructure()
                v = sparse(ones(nam - 1, nam));
            end                
            
            output.convertFromState  = @(state) convertFromState(state);
            output.convertToState    = @(x) convertToState(x);

            funcs.objective = @(x) computeVoltage(x);
            funcs.gradient  = @(x) computeDerVoltage(x);
            
            funcs.constraints       = @(x) computeConstraints(x);
            funcs.jacobian          = @(x) computeDerConstraints(x);
            funcs.jacobianstructure = @() getConstraintJacobianStructure;

            lb = [];
            ub = [];
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                lb = vertcat(lb, ...
                             model.(elde).thetaMins);
                ub = vertcat(ub, ...
                             model.(elde).thetaMaxs);
            end

            cl = zeros(nam - 1, 1);
            cu = cl;

            options = struct('lb', lb, ...
                             'ub', ub, ...
                             'cl', cl, ...
                             'cu', cu);

        end

        function [stateStart, stateEnd] = computeExtremalStates(ecs, stateInit, varargin)

            opt = struct('verbose', false);
            opt = merge_options(opt, varargin{:});
            
            [funcs, options, extras] = ecs.setupIpOpt('max');
            options.ipopt.hessian_approximation = 'limited-memory';
            if ~opt.verbose
                options.ipopt.print_level = 0;
            end
            x0 = extras.convertFromState(stateInit);
            [x, info] = ipopt(x0, funcs, options);
            stateStart = extras.convertToState(x);
            
            [funcs, options, extras] = ecs.setupIpOpt('min');
            options.ipopt.hessian_approximation = 'limited-memory';
            if ~opt.verbose
                options.ipopt.print_level = 0;
            end
            x0 = extras.convertFromState(stateInit);
            [x, info] = ipopt(x0, funcs, options);
            stateEnd = extras.convertToState(x);

        end
        
    end

end
