classdef OxideMembraneCell < BaseModel
    
    properties
        
        % Temperature
        T
        % Structure with physical constants
        constants

        Anode
        Cathode
        Electrolyte
        Control
        
        couplingTerms
        couplingnames

        dx
        farea

        primaryVarNames
        equationVarNames
        residualFuncCallList
        funcCallList        
        
    end
    
    methods
        
        function model = OxideMembraneCell(paramobj)

            model = model@BaseModel();

            fdnames = {'T'    , ...
                       'dx'   , ...
                       'farea', ...
                       'couplingTerms'};
            model = dispatchParams(model, paramobj, fdnames);

            model.Anode       = OxideMembraneElectrode(paramobj.Anode);
            model.Cathode     = OxideMembraneElectrode(paramobj.Cathode);
            model.Electrolyte = OxideMembraneElectrolyte(paramobj.Electrolyte);
            model.Control     = OxideMembraneControl(paramobj.Control);
            
            % setup couplingNames
            model.couplingnames = cellfun(@(x) x.name, model.couplingTerms, 'uniformoutput', false);

            % setup standard physical constants
            model.constants = PhysicalConstants();

        end
        
        function model = registerVarAndPropfuncNames(model)
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            fn = @OxideMembraneCell.setupO2Sources;
            inputnames = {{an, 'jO2m'}, {ct, 'jO2m'}};
            model = model.registerPropFunction({{elyte, 'sourceO2'}, fn, inputnames});
            
            fn = @OxideMembraneCell.setupElSources;
            inputnames = {{an, 'jEl'}, {ct, 'jEl'}};
            model = model.registerPropFunction({{elyte, 'sourceEl'}, fn, inputnames});

            fn = @OxideMembraneCell.updateAnodeJO2Equation;
            inputnames = {{elyte, 'phi'}, {an, 'phi'}, {an, 'jO2m'}};
            model = model.registerPropFunction({{an, 'jO2mEquation'}, fn, inputnames});

            fn = @OxideMembraneCell.updateAnodeJElEquation;
            inputnames = {{elyte, 'gradPhiCoef'}              , ...
                          VarName({elyte}, 'gradConcCoefs', 2), ...
                          {elyte, 'phi'}                      , ...
                          {elyte, 'ce'}                       , ...
                          {elyte, 'ch'}                       , ...
                          {an, 'phi'}                         , ...
                          {an, 'ce'}                          , ...
                          {an, 'ch'}                          , ...
                          {an, 'jEl'}};
            model = model.registerPropFunction({{an, 'jElEquation'}, fn, inputnames});

            fn = @OxideMembraneCell.updateCathodeJO2Equation;
            inputnames = {{elyte, 'phi'}, {ct, 'phi'}, {ct, 'jO2m'}};
            model = model.registerPropFunction({{ct, 'jO2mEquation'}, fn, inputnames});

            fn = @OxideMembraneCell.updateCathodeJElEquation;
            inputnames = {{elyte, 'gradPhiCoef'}              , ...
                          VarName({elyte}, 'gradConcCoefs', 2), ...
                          {elyte, 'ce'}                       , ...
                          {elyte, 'ch'}                       , ...
                          {ct, 'phi'}                         , ...
                          {ct, 'ce'}                          , ...
                          {ct, 'ch'}                          , ...
                          {ct, 'jEl'}};
            model = model.registerPropFunction({{ct, 'jElEquation'}, fn, inputnames});
            
            fn = @OxideMembraneCell.updateFromControl;
            inputnames = {{ctrl, 'U'}, {ctrl, 'I'}};
            model = model.registerPropFunction({{an, 'pi'}, fn, inputnames});                        
            model = model.registerPropFunction({{an, 'jO2m'}, fn, inputnames});

            fn = @OxideMembraneCell.setupCathodeBoundary;
            inputnames = {};
            model = model.registerPropFunction({{ct, 'phi'}, fn, inputnames});
                        
            inputnames = {};
            fn = @OxideMembraneCell.updateControl;
            fn = {fn, @(propfunction) PropFunction.drivingForceFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({{ctrl, 'ctrlVal'}, fn, inputnames});            
            model = model.registerPropFunction({{elyte, 'alpha'}, fn, inputnames});
            
        end

        function state = updateControl(model, state, drivingForces)
            
            ctrl = "Control";
            elyte = 'Electrolyte';
            
            time = state.time;
            [ctrlVal, alpha] = drivingForces.src(time);
            
            state.(ctrl).ctrlVal = ctrlVal;
            state.(elyte).alpha  = alpha;
            
        end
        
        function state = setupO2Sources(model, state)

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';

            coupterms = model.couplingTerms;
            coupnames = model.couplingnames;
                        
            jO2mAnode   = state.(an).jO2m;
            jO2mCathode = state.(ct).jO2m;
            
            sourceO2 = 0*state.(elyte).phi; % initialize AD for sourceO2
            
            % Anode part
            coupterm = getCoupTerm(coupterms, 'Anode-Electrolyte', coupnames);
            ccs = coupterm.couplingcells(:, 2);
            sourceO2(ccs) = jO2mAnode;
            
            % Anode part
            coupterm = getCoupTerm(coupterms, 'Cathode-Electrolyte', coupnames);
            ccs = coupterm.couplingcells(:, 2);
            sourceO2(ccs) = jO2mCathode;

            state.(elyte).sourceO2 = sourceO2;
            
        end

        function state = setupElSources(model, state)
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';

            coupterms = model.couplingTerms;
            coupnames = model.couplingnames;
            
            jElAnode   = state.(an).jEl;
            jElCathode = state.(ct).jEl;
            
            sourceEl = 0*state.(elyte).phi; % initialize AD for sourceEl
            
            % Anode part
            coupterm = getCoupTerm(coupterms, 'Anode-Electrolyte', coupnames);
            ccs = coupterm.couplingcells(:, 2);
            sourceEl(ccs) = jElAnode;
            
            % Anode part
            coupterm = getCoupTerm(coupterms, 'Cathode-Electrolyte', coupnames);
            ccs = coupterm.couplingcells(:, 2);
            sourceEl(ccs) = jElCathode;

            state.(elyte).sourceEl = sourceEl;
            
        end

        function state = updateFromControl(model, state)

            an   = 'Anode';
            ctrl = 'Control';
            
            state.(an).pi   = state.(ctrl).U;
            state.(an).jO2m = state.(ctrl).I;
            
        end

        function state = setupCathodeBoundary(model, state)

            ct   = 'Cathode';

            state.(ct).phi = 0;
            
        end

        function state = updateCathodeJElEquation(model, state)

            state = model.updateJElEquation(state, 'Cathode');
            
        end

        function state = updateAnodeJElEquation(model, state)

            state = model.updateJElEquation(state, 'Anode');
            
        end

        function state = updateCathodeJO2Equation(model, state)

            state = model.updateJO2Equation(state, 'Cathode');
            
        end

        function state = updateAnodeJO2Equation(model, state)

            state = model.updateJO2Equation(state, 'Anode');
            
        end
        
        function state = updateJElEquation(model, state, elde)

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';

            op    = model.(elyte).operators;
            cinds = model.(elyte).compinds;

            coupterms = model.couplingTerms;
            coupnames = model.couplingnames;

            switch elde
              case an
                coupTerm = getCoupTerm(coupterms, 'Anode-Electrolyte', coupnames);
              case ct
                coupTerm = getCoupTerm(coupterms, 'Cathode-Electrolyte', coupnames);
              otherwise
                error('Coupling term not found');
            end

            ccs = coupTerm.couplingcells;
            cfs = coupTerm.couplingfaces;

            phiElde  = state.(elde).phi;
            ceElde   = state.(elde).ce;
            chElde   = state.(elde).ch;
            jEl      = state.(elde).jEl;
            phiElyte = state.(elyte).phi;
            ceElyte  = state.(elyte).ce;
            chElyte  = state.(elyte).ch;
            gPhiC    = state.(elyte).gradPhiCoef;
            gConcCs  = state.(elyte).gradConcCoefs;
            
            t    = op.harmFaceBC(gConcCs{cinds.ch}, cfs(:, 2));
            jch  = t.*(chElde(ccs(:, 1)) - chElyte(ccs(:, 2)));

            t    = op.harmFaceBC(gConcCs{cinds.ce}, cfs(:, 2));
            jce  = t.*(ceElde(ccs(:, 1)) - ceElyte(ccs(:, 2)));

            t    = op.harmFaceBC(gPhiC, cfs(:, 2));
            jphi = t.*(phiElde(ccs(:, 1)) - phiElyte(ccs(:, 2)));

            state.(elde).jElEquation = jEl - (jch - jce + jphi);

        end

        function state = updateJO2Equation(model, state, elde)

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';

            op = model.(elyte).operators;

            coupterms = model.couplingTerms;
            coupnames = model.couplingnames;
            
            switch elde
              case an
                coupTerm = getCoupTerm(coupterms, 'Anode-Electrolyte', coupnames);
              case ct
                coupTerm = getCoupTerm(coupterms, 'Cathode-Electrolyte', coupnames);
              otherwise
                error('Coupling term not found');
            end

            ccs = coupTerm.couplingcells;
            cfs = coupTerm.couplingfaces;
            
            sigmaO2  = model.(elyte).sigmaO2;

            phiElyte = state.(elyte).phi;
            phiElde  = state.(elde).phi;
            jO2m     = state.(elde).jO2m;
            
            tcoef = op.halfTransBC(cfs(:, 2));
            
            % Note plus-sign between jO2m and sigmaO2, because O2m is negatively charged.
            state.(elde).jO2mEquation = jO2m + sigmaO2*tcoef.*(phiElde(ccs(:, 1)) - phiElyte(ccs(:, 2)));

        end
        

        function initState = setupInitialState(model)
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';
            
            T = model.T;
            c = model.constants;
            
            %% Anode initialization
            % jEl is setup afterwards
            
            Eocp  = model.(an).Eocp;
            mu0   = model.(an).muEl0;
            K     = model.(an).Keh;
            cinds = model.(an).compinds;
            
            initState.(an).jO2m            = 0;
            initState.(an).j               = 0;
            initState.(an).logcs{cinds.ce} = - (c.F*Eocp + mu0)/(c.R*T);
            initState.(an).logcs{cinds.ch} = log(K) + (c.F*Eocp + mu0)/(c.R*T); % sum of log should be equal to log(K)
            initState.(an).phi             = 0;
            
            %% Cathode initialization
            %  jEl is setup afterwards
            
            Eocp  = model.(ct).Eocp;
            mu0   = model.(ct).muEl0;
            K     = model.(ct).Keh;
            cinds = model.(ct).compinds;

            initState.(ct).jO2m            = 0;
            initState.(ct).j               = 0;
            initState.(ct).logcs{cinds.ce} = -(c.F*Eocp + mu0)/(c.R*T);
            initState.(ct).logcs{cinds.ch} = log(K) + (c.F*Eocp + mu0)/(c.R*T); % sum of log should be equal to log(K)
            initState.(ct).pi              = model.(ct).Eocp;

            %% Electrolyte initialization
            %
            
            nc = model.(elyte).G.cells.num;
            initState.(elyte).phi = zeros(nc, 1);

            % we initiate with a linear interpolation

            lcan = initState.(an).logcs{cinds.ce};
            lcct = initState.(ct).logcs{cinds.ce};

            G = model.(elyte).G;

            cc   = G.cells.centroids(:, 1);
            ndan = G.nodes.coords(1, 1);
            ndct = G.nodes.coords(end, 1);

            lc = interp1([ndan; ndct], [lcan; lcct], cc);
            
            initState.(elyte).logcs{cinds.ce} = lc;
            initState.(elyte).logcs{cinds.ch} = log(K) - lc; % sum of log should be equal to log(K)

            initState.time = 0;

            %% Initiate jEl at electrodes using jElEquation
            %
            
            function [ctrlVal, alpha] = src(time)
                ctrlVal = 0;
                alpha = 0;
            end

            drivingForces.src = @(time) src(time);

            initState.(an).jEl = 0; 
            initState = model.evalVarName(initState, {an, 'jElEquation'}, {{'drivingForces', drivingForces}});
            initState.(an).jEl = -initState.(an).jElEquation;
            
            initState.(ct).jEl = 0;
            initState = model.evalVarName(initState, {ct, 'jElEquation'}, {{'drivingForces', drivingForces}});
            initState.(ct).jEl = -initState.(ct).jElEquation;

            %% Initiate Control
            %

            initState.(ctrl).I = 0;
            initState.(ctrl).U = model.(an).Eocp;

        end

        function state = addVariables(model, state, drivingForces)
                        
            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end
            
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)

            opts = struct('ResOnly', false, 'iteration', 0, 'reverseMode', false); 
            opts = merge_options(opts, varargin{:});
            
            state.time = state0.time + dt;
            
            if(not(opts.ResOnly) && not(opts.reverseMode))
                state = model.initStateAD(state);
            elseif(opts.reverseMode)
                disp('No AD initatlization in equation old style')
                state0 = model.initStateAD(state0);
            else
                assert(opts.ResOnly);
            end
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';

            %% We call the assembly equations ordered from the graph
            
            funcCallList = model.residualFuncCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end

            eqnvarnames = model.equationVarNames;
            neq = numel(eqnvarnames);

            for ieq = 1 : neq

                varname = eqnvarnames{ieq};

                eqs{ieq}   = model.getProp(state, varname);
                names{ieq} = strjoin(varname, '.');

            end

            types = repmat({'cells'}, 1, numel(names));
            
            primaryVars = model.getPrimaryVariables();
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end
        
        function primaryvarnames = getPrimaryVariableNames(model)

            primaryvarnames = model.primaryVarNames;
            
        end
        
        function forces = getValidDrivingForces(model)
            
            forces = getValidDrivingForces@PhysicalModel(model);
            forces.src = [];
            forces.alpha = [];
            
        end

        function model = validateModel(model, varargin)

            if isempty(model.computationalGraph)
                model = model.setupComputationalGraph();
            end

            cgit = model.computationalGraph;
            
            model.primaryVarNames      = cgt.getPrimaryVariableNames();
            model.equationVarNames     = cgt.getEquationVariableNames();
            model.residualFuncCallList = cgt.getOrderedFunctionCallList();
            model.funcCallList         = cgt.getOrderedFunctionCallList('removeExtraVariables', false);
            
        end
        
    end
    
end
