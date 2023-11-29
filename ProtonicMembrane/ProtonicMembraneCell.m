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

        standalone

    end

    methods

        function model = ProtonicMembraneCell(paramobj)

            model = model@BaseModel();

            fdnames = {'T'       , ...
                       'couplingTerms'};

            model = dispatchParams(model, paramobj, fdnames);

            model.Anode       = ProtonicMembraneAnode(paramobj.Anode);
            model.Cathode     = ProtonicMembraneCathode(paramobj.Cathode);
            model.Electrolyte = ProtonicMembraneElectrolyte(paramobj.Electrolyte);
            model.Control     = ProtonicMembraneControl(paramobj.Control);


            cps = model.couplingTerms;
            assert(all(size(cps{1}.couplingcells) == size(cps{2}.couplingcells)), ...
                   'The coupling terms appear to be badly setup');
            model.Anode.N   = size(cps{1}.couplingcells, 1);
            model.Cathode.N = size(cps{1}.couplingcells, 1);

            % setup couplingNames
            model.couplingnames = cellfun(@(x) x.name, model.couplingTerms, 'uniformoutput', false);

            % setup standard physical constants
            model.constants = PhysicalConstants();

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';

            sigmaHp = model.(elyte).sigma_prot;
            sigmaEl = model.(elyte).sigma_n0;
            phi0    = abs(model.(an).E_0 - model.(ct).E_0); % characteristic voltage
            T       = model.(elyte).operators.T_all(1);

            sHp = T*sigmaHp*phi0;
            sEl = T*sigmaEl*phi0;

            model.scalings =  {{{elyte, 'massConsHp'}     , sHp}, ...
                               {{elyte, 'chargeConsEl'}   , sEl}, ...
                               {{an   , 'chargeCons'}     , sEl}, ...
                               {{an   , 'jElEquation'}    , sEl}, ...
                               {{an   , 'jHpEquation'}    , sHp}, ...
                               {{ct   , 'chargeCons'}     , sEl}, ...
                               {{ct   , 'jElEquation'}    , sEl}, ...
                               {{ct   , 'jHpEquation'}    , sHp}, ...
                               {{ctrl , 'controlEquation'}, sEl}, ...
                               {{'anodeChargeCons'}       , sEl}};

        end

        function model = registerVarAndPropfuncNames(model)

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {};
            varnames{end + 1} = 'anodeChargeCons';

            model = model.registerVarNames(varnames);

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
            inputnames = {{ctrl, 'U'}};
            model = model.registerPropFunction({{an, 'pi'}, fn, inputnames});

            fn = @ProtonicMembraneCell.updateAnodeChargeCons;
            inputnames = {{ctrl, 'I'}, {an, 'j'}};
            model = model.registerPropFunction({'anodeChargeCons', fn, inputnames});

            fn = @ProtonicMembraneCell.setupCathodeBoundary;
            inputnames = {};
            model = model.registerPropFunction({{ct, 'phi'}, fn, inputnames});

            inputnames = {};
            fn = @ProtonicMembraneCell.updateControl;
            fn = {fn, @(propfunction) PropFunction.drivingForceFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({{ctrl, 'ctrlVal'}, fn, inputnames});
            model = model.registerPropFunction({{elyte, 'alpha'}, fn, inputnames});

        end

        function state = updateAnodeChargeCons(model, state)

            an   = 'Anode';
            ctrl = 'Control';

            state.anodeChargeCons = sum(state.(an).j) - state.(ctrl).I;

        end

        function state = updateControl(model, state, drivingForces)

            ctrl  = "Control";
            an    = 'Anode';
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

            N = model.(an).N;

            onevec = ones(N, 1);
            state.(an).pi = state.(ctrl).U.*onevec;

        end

        function state = setupCathodeBoundary(model, state)

            ct   = 'Cathode';

            state.(ct).phi = 0*ones(model.(ct).N, 1);

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

            state.(elde).jElEquation = jEl - T.*(piElde(ccs(:, 1)) - piElyte(ccs(:, 2)));

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

            state.(elde).jHpEquation = jHp - T.*(phiElde(ccs(:, 1)) - phiElyte(ccs(:, 2)));

        end


        function initState = setupInitialState(model)

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';

            onevec = ones(model.(an).N, 1);

            initState.(an).phi = 0*onevec;
            initState.(an).jHp = 0*onevec;
            initState.(an).jEl = 0*onevec;
            initState.(an).j   = 0*onevec;

            nc = model.(elyte).G.cells.num;
            initState.(elyte).pi  = zeros(nc, 1);
            initState.(elyte).phi = zeros(nc, 1);

            onevec = ones(model.(ct).N, 1);

            initState.(ct).pi  = 0*onevec;
            initState.(ct).jEl = 0*onevec;
            initState.(ct).jHp = 0*onevec;
            initState.(ct).j   = 0*onevec;

            initState.(ctrl).I = 0;

            initState = model.evalVarName(initState, {an, 'Eocv'});
            initState.(ctrl).U = initState.(an).Eocv(1);

            initState.time = 0;

        end

        function model = setupStandAlone(model)

            model.standalone = true;

            shortNames = {'Electrolyte', 'elyte';
                          'Anode'      , 'an';
                          'Cathode'    , 'ct';
                          'Control'    , 'ctrl'};

            model = model.equipModelForComputation('shortNames', shortNames);

        end


        function state = addVariables(model, state, drivingForces)

            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end

        end

        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@BaseModel(model);
            forces.src = [];
            forces.alpha = [];

        end

    end

end
