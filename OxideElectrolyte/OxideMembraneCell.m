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

            model.Anode       = OxideMembraneAnode(paramobj.Anode);
            model.Cathode     = OxideMembraneCathode(paramobj.Cathode);
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
            inputnames = {{an, 'jO2'}, {ct, 'jO2'}};
            model = model.registerPropFunction({{elyte, 'sourceO2'}, fn, inputnames});
            
            fn = @OxideMembraneCell.setupElSources;
            inputnames = {{an, 'jEl'}, {ct, 'jEl'}};
            model = model.registerPropFunction({{elyte, 'sourceEl'}, fn, inputnames});

            fn = @OxideMembraneCell.updateAnodeJO2Equation;
            inputnames = {{elyte, 'phi'}, {an, 'phi'}, {elyte, 'sigmaO2'}, {an, 'jO2'}};
            model = model.registerPropFunction({{an, 'jO2Equation'}, fn, inputnames});

            fn = @OxideMembraneCell.updateAnodeJElEquation;
            inputnames = {{elyte, 'phi'}, {an, 'phi'}, {elyte, 'sigmaEl'}, {an, 'jEl'}};
            model = model.registerPropFunction({{an, 'jElEquation'}, fn, inputnames});

            fn = @OxideMembraneCell.updateCathodeJO2Equation;
            inputnames = {{elyte, 'phi'}, {ct, 'phi'}, {elyte, 'sigmaO2'}, {ct, 'jO2'}};
            model = model.registerPropFunction({{ct, 'jO2Equation'}, fn, inputnames});

            fn = @OxideMembraneCell.updateCathodeJElEquation;
            inputnames = {{elyte, 'pi'}, {ct, 'pi'}, {elyte, 'sigmaEl'}, {ct, 'jEl'}};
            model = model.registerPropFunction({{ct, 'jElEquation'}, fn, inputnames});
            
            fn = @OxideMembraneCell.updateFromControl;
            inputnames = {{ctrl, 'U'}, {ctrl, 'I'}};
            model = model.registerPropFunction({{an, 'pi'}, fn, inputnames});                        
            model = model.registerPropFunction({{an, 'j'}, fn, inputnames});

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
                        
            jO2Anode   = state.(an).jO2;
            jO2Cathode = state.(ct).jO2;
            
            sourceO2 = 0*state.(elyte).phi; % initialize AD for sourceO2
            
            % Anode part
            coupterm = getCoupTerm(coupterms, 'Anode-Electrolyte', coupnames);
            ccs = coupterm.couplingcells(:, 2);
            sourceO2(ccs) = jO2Anode;
            
            % Anode part
            coupterm = getCoupTerm(coupterms, 'Cathode-Electrolyte', coupnames);
            ccs = coupterm.couplingcells(:, 2);
            sourceO2(ccs) = jO2Cathode;

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
            
            state.(an).pi  = state.(ctrl).U; 
            state.(an).j   = state.(ctrl).I;
            
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
            
            sigmaEl = state.(elyte).sigmaEl;
            piElyte = state.(elyte).pi;
            piElde  = state.(elde).pi;
            jEl     = state.(elde).jEl;

            T = op.harmFaceBC(sigmaEl, cfs(:, 2));

            state.(elde).jElEquation = jEl - T*(piElde(ccs(:, 1)) - piElyte(ccs(:, 2)));

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
            
            sigmaO2  = state.(elyte).sigmaO2;
            phiElyte = state.(elyte).phi;
            phiElde  = state.(elde).phi;
            jO2      = state.(elde).jO2;

            T = op.harmFaceBC(sigmaO2, cfs(:, 2));

            state.(elde).jO2Equation = jO2 - T*(phiElde(ccs(:, 1)) - phiElyte(ccs(:, 2)));

        end
        

        function initState = setupInitialState(model)
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';
            
            
            initState.(an).pi  = model.(an).Eocp;
            initState.(an).phi = 0;
            initState.(an).jO2 = 0;
            initState.(an).jEl = 0;
            
            nc = model.(elyte).G.cells.num;
            initState.(elyte).pi  = zeros(nc, 1);
            initState.(elyte).phi = zeros(nc, 1);

            initState.(ct).pi  = 0;
            initState.(ct).jEl = 0;
            initState.(ct).jO2 = 0;
            initState.(ct).j   = 0;

            initState.(ctrl).I = 0;
            initState.(ctrl).U = model.(an).Eocp;

            initState.time = 0;
            
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

            %% We call the assembly equations ordered from the graph
            
            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';

            eqs = {};

            dx      = model.dx;
            farea   = model.farea;
            sigmaO2 = model.(elyte).sigmaO2;
            sigmaEl = model.(elyte).sigmaN_0;
            phi0    = 1;
            
            sO2 = 1/(farea*sigmaO2*phi0/dx);
            sEl = 1/(farea*sigmaEl*phi0/dx);
            
            eqs{end + 1} = sO2*state.(elyte).massConsO2;
            eqs{end + 1} = sEl*state.(elyte).chargeConsEl;
            eqs{end + 1} = sEl*state.(an).chargeCons;
            eqs{end + 1} = sEl*state.(an).jElEquation;
            eqs{end + 1} = sO2*state.(an).jO2Equation;
            eqs{end + 1} = sEl*state.(ct).chargeCons;
            eqs{end + 1} = sEl*state.(ct).jElEquation;
            eqs{end + 1} = sO2*state.(ct).jO2Equation;
            eqs{end + 1} = sEl*state.(ctrl).controlEquation;
            
            names = {'elyte_massConsO2'  , ...
                     'elyte_chargeConsEl', ...
                     'an_chargeCons'     , ...
                     'an_jElEquation'    , ...
                     'an_jO2Equation'    , ...
                     'ct_chargeCons'     , ...
                     'ct_jElEquation'    , ...
                     'ct_jO2Equation'    , ...
                     'ctrl_controlEquation'};

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

            cgt = model.computationalGraph;
            
            model.primaryVarNames = cgt.getPrimaryVariableNames();
            model.funcCallList    = cgt.getOrderedFunctionCallList();
            
        end
        
    end
    
end
