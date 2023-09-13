classdef ProtonicMembraneCell < BaseModel
    
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

        primaryVarNames
        funcCallList
        
    end
    
    methods
        
        function model = ProtonicMembraneCell(paramobj)

            model = model@BaseModel();

            fdnames = {'T', ...
                       'couplingTerms'};
            model = dispatchParams(model, paramobj, fdnames);

            model.Anode       = ProtonicMembraneAnode(paramobj.Anode);
            model.Cathode     = ProtonicMembraneCathode(paramobj.Cathode);
            model.Electrolyte = ProtonicMembraneElectrolyte(paramobj.Electrolyte);
            model.Control     = ProtonicMembraneControl(paramobj.Control);
            
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
            
            fn = @ProtonicMembraneCell.setupHpSources;
            inputnames = {{an, 'jHp'}, {ct, 'jHp'}};
            model = model.registerPropFunction({{elyte, 'sourceHp'}, fn, inputnames});
            
            fn = @ProtonicMembraneCell.setupElSources;
            inputnames = {{an, 'jEl'}, {ct, 'jEl'}};
            model = model.registerPropFunction({{elyte, 'sourceEl'}, fn, inputnames});

            fn = @ProtonicMembraneCell.updateAnodeJHpEquation;
            inputnames = {{elyte, 'phi'}, {an, 'phi'}, {elyte, 'sigmaHp'}, {an, 'jHp'}};
            model = model.registerPropFunction({{an, 'jHpEquation'}, fn, inputnames});

            fn = @ProtonicMembraneCell.updateAnodeJElEquation;
            inputnames = {{elyte, 'phi'}, {an, 'phi'}, {elyte, 'sigmaEl'}, {an, 'jEl'}};
            model = model.registerPropFunction({{an, 'jElEquation'}, fn, inputnames});

            fn = @ProtonicMembraneCell.updateCathodeJHpEquation;
            inputnames = {{elyte, 'phi'}, {ct, 'phi'}, {elyte, 'sigmaHp'}, {ct, 'jHp'}};
            model = model.registerPropFunction({{ct, 'jHpEquation'}, fn, inputnames});

            fn = @ProtonicMembraneCell.updateCathodeJElEquation;
            inputnames = {{elyte, 'pi'}, {ct, 'pi'}, {elyte, 'sigmaEl'}, {ct, 'jEl'}};
            model = model.registerPropFunction({{ct, 'jElEquation'}, fn, inputnames});
            
            fn = @ProtonicMembraneCell.updateFromControl;
            inputnames = {{ctrl, 'U'}, {ctrl, 'I'}};
            model = model.registerPropFunction({{an, 'pi'}, fn, inputnames});                        
            model = model.registerPropFunction({{an, 'j'}, fn, inputnames});

            fn = @ProtonicMembraneCell.setupCathodeBoundary;
            inputnames = {};
            model = model.registerPropFunction({{ct, 'phi'}, fn, inputnames});
                        
            inputnames = {};
            fn = @ProtonicMembraneCell.updateControl;
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
        
        function state = setupHpSources(model, state)

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';

            coupterms = model.couplingTerms;
            coupnames = model.couplingnames;
                        
            jHpAnode   = state.(an).jHp;
            jHpCathode = state.(ct).jHp;
            
            sourceHp = 0*state.(elyte).phi; % initialize AD for sourceHp
            
            % Anode part
            coupterm = getCoupTerm(coupterms, 'Anode-Electrolyte', coupnames);
            ccs = coupterm.couplingcells(:, 2);
            sourceHp(ccs) = jHpAnode;
            
            % Anode part
            coupterm = getCoupTerm(coupterms, 'Cathode-Electrolyte', coupnames);
            ccs = coupterm.couplingcells(:, 2);
            sourceHp(ccs) = jHpCathode;

            state.(elyte).sourceHp = sourceHp;
            
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

        function state = updateCathodeJHpEquation(model, state)

            state = model.updateJHpEquation(state, 'Cathode');
            
        end

        function state = updateAnodeJHpEquation(model, state)

            state = model.updateJHpEquation(state, 'Anode');
            
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

        function state = updateJHpEquation(model, state, elde)

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
            
            sigmaHp  = state.(elyte).sigmaHp;
            phiElyte = state.(elyte).phi;
            phiElde  = state.(elde).phi;
            jHp      = state.(elde).jHp;

            T = op.harmFaceBC(sigmaHp, cfs(:, 2));

            state.(elde).jHpEquation = jHp - T*(phiElde(ccs(:, 1)) - phiElyte(ccs(:, 2)));

        end
        

        function initState = setupInitialState(model)
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';
            
            
            initState.(an).pi  = model.(an).Eocp;
            initState.(an).phi = 0;
            initState.(an).jHp = 0;
            initState.(an).jEl = 0;
            
            nc = model.(elyte).G.cells.num;
            initState.(elyte).pi  = zeros(nc, 1);
            initState.(elyte).phi = zeros(nc, 1);

            initState.(ct).pi  = 0;
            initState.(ct).jEl = 0;
            initState.(ct).jHp = 0;
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
            eqs{end + 1} = state.(elyte).massConsHp;
            eqs{end + 1} = state.(elyte).chargeConsEl;
            eqs{end + 1} = state.(an).chargeCons;
            eqs{end + 1} = state.(an).jElEquation;
            eqs{end + 1} = state.(an).jHpEquation;
            eqs{end + 1} = state.(ct).chargeCons;
            eqs{end + 1} = state.(ct).jElEquation;
            eqs{end + 1} = state.(ct).jHpEquation;
            eqs{end + 1} = state.(ctrl).controlEquation;
            
            names = {'elyte_massConsHp'  , ...
                     'elyte_chargeConsEl', ...
                     'an_chargeCons'     , ...
                     'an_jElEquation'    , ...
                     'an_jHpEquation'    , ...
                     'ct_chargeCons'     , ...
                     'ct_jElEquation'    , ...
                     'ct_jHpEquation'    , ...
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
