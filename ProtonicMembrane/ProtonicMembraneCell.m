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

    end

    methods

        function model = ProtonicMembraneCell(inputparams)

            model = model@BaseModel();

            fdnames = {'G', ...
                       'T', ...
                       'couplingTerms'};

            model = dispatchParams(model, inputparams, fdnames);

            model.Anode       = ProtonicMembraneAnode(inputparams.Anode);
            model.Cathode     = ProtonicMembraneCathode(inputparams.Cathode);
            model.Electrolyte = ProtonicMembraneElectrolyte(inputparams.Electrolyte);
            model.Control     = ProtonicMembraneControl(inputparams.Control);


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
            T       = model.(elyte).G.getTrans();
            T       = mean(T);

            sHp = T*sigmaHp*phi0;
            sEl = T*sigmaEl*phi0;

            model.scalings =  {{{elyte, 'massConsHp'}     , sHp}, ...
                               {{elyte, 'chargeConsEl'}   , sEl}, ...
                               {{an   , 'chargeCons'}     , sEl}, ...
                               {{an   , 'iElEquation'}    , sEl}, ...
                               {{an   , 'iHpEquation'}    , sHp}, ...
                               {{ct   , 'chargeCons'}     , sEl}, ...
                               {{ct   , 'iElEquation'}    , sEl}, ...
                               {{ct   , 'iHpEquation'}    , sHp}, ...
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
            varnames{end + 1} = 'time';
            
            model = model.registerVarNames(varnames);

            fn = @ProtonicMembraneCell.setupHpSources;
            inputnames = {{an, 'iHp'}, {ct, 'iHp'}};
            model = model.registerPropFunction({{elyte, 'sourceHp'}, fn, inputnames});

            fn = @ProtonicMembraneCell.setupElSources;
            inputnames = {{an, 'iEl'}, {ct, 'iEl'}};
            model = model.registerPropFunction({{elyte, 'sourceEl'}, fn, inputnames});

            fn = @ProtonicMembraneCell.updateAnodeIHpEquation;
            inputnames = {{elyte, 'phi'}, {an, 'phi'}, {elyte, 'sigmaHp'}, {an, 'iHp'}};
            model = model.registerPropFunction({{an, 'iHpEquation'}, fn, inputnames});

            fn = @ProtonicMembraneCell.updateAnodeIElEquation;
            inputnames = {{elyte, 'pi'}, {an, 'pi'}, {elyte, 'sigmaEl'}, {an, 'iEl'}};
            model = model.registerPropFunction({{an, 'iElEquation'}, fn, inputnames});

            fn = @ProtonicMembraneCell.updateCathodeIHpEquation;
            inputnames = {{elyte, 'phi'}, {ct, 'phi'}, {elyte, 'sigmaHp'}, {ct, 'iHp'}};
            model = model.registerPropFunction({{ct, 'iHpEquation'}, fn, inputnames});

            fn = @ProtonicMembraneCell.updateCathodeIElEquation;
            inputnames = {{elyte, 'pi'}, {ct, 'pi'}, {elyte, 'sigmaEl'}, {ct, 'iEl'}};
            model = model.registerPropFunction({{ct, 'iElEquation'}, fn, inputnames});

            fn = @ProtonicMembraneCell.updateFromControl;
            inputnames = {{ctrl, 'U'}};
            model = model.registerPropFunction({{an, 'pi'}, fn, inputnames});

            fn = @ProtonicMembraneCell.updateAnodeChargeCons;
            inputnames = {{ctrl, 'I'}, {an, 'i'}};
            model = model.registerPropFunction({'anodeChargeCons', fn, inputnames});

            fn = @ProtonicMembraneCell.setupCathodeBoundary;
            inputnames = {};
            model = model.registerPropFunction({{ct, 'phi'}, fn, inputnames});

            inputnames = {'time'};
            fn = @ProtonicMembraneCell.updateControl;
            fn = {fn, @(propfunction) PropFunction.drivingForceFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({{ctrl, 'ctrlVal'}, fn, inputnames});
            model = model.registerPropFunction({{elyte, 'alpha'}, fn, inputnames});

            model = model.setAsStaticVarName('time');
            
        end

        function state = updateAnodeChargeCons(model, state)

            an   = 'Anode';
            ctrl = 'Control';

            state.anodeChargeCons = sum(state.(an).i) - state.(ctrl).I;

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

            iHpAnode   = state.(an).iHp;
            iHpCathode = state.(ct).iHp;

            sourceHp = 0*state.(elyte).phi; % initialize AD for sourceHp

            % Anode part
            coupterm = getCoupTerm(coupterms, 'Anode-Electrolyte', coupnames);
            ccs = coupterm.couplingcells(:, 2);
            sourceHp(ccs) = iHpAnode;

            % Anode part
            coupterm = getCoupTerm(coupterms, 'Cathode-Electrolyte', coupnames);
            ccs = coupterm.couplingcells(:, 2);
            sourceHp(ccs) = iHpCathode;

            state.(elyte).sourceHp = sourceHp;

        end

        function state = setupElSources(model, state)

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';

            coupterms = model.couplingTerms;
            coupnames = model.couplingnames;

            iElAnode   = state.(an).iEl;
            iElCathode = state.(ct).iEl;

            sourceEl = 0*state.(elyte).phi; % initialize AD for sourceEl

            % Anode part
            coupterm = getCoupTerm(coupterms, 'Anode-Electrolyte', coupnames);
            ccs = coupterm.couplingcells(:, 2);
            sourceEl(ccs) = iElAnode;

            % Anode part
            coupterm = getCoupTerm(coupterms, 'Cathode-Electrolyte', coupnames);
            ccs = coupterm.couplingcells(:, 2);
            sourceEl(ccs) = iElCathode;

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

        function state = updateCathodeIElEquation(model, state)

            state = model.updateIElEquation(state, 'Cathode');

        end

        function state = updateAnodeIElEquation(model, state)

            state = model.updateIElEquation(state, 'Anode');

        end

        function state = updateCathodeIHpEquation(model, state)

            state = model.updateIHpEquation(state, 'Cathode');

        end

        function state = updateAnodeIHpEquation(model, state)

            state = model.updateIHpEquation(state, 'Anode');

        end


        function state = updateIElEquation(model, state, elde)

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
            iEl     = state.(elde).iEl;

            T = model.(elyte).G.getTransBcHarmFace(sigmaEl, cfs(:, 2));

            state.(elde).iElEquation = iEl - T.*(piElde(ccs(:, 1)) - piElyte(ccs(:, 2)));

        end

        function state = updateIHpEquation(model, state, elde)

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
            iHp      = state.(elde).iHp;

            T = model.(elyte).G.getTransBcHarmFace(sigmaHp, cfs(:, 2));
            
            state.(elde).iHpEquation = iHp - T.*(phiElde(ccs(:, 1)) - phiElyte(ccs(:, 2)));

        end


        function initState = setupInitialState(model)

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';

            onevec = ones(model.(an).N, 1);

            initState.(an).phi = 0*onevec;
            initState.(an).iHp = 0*onevec;
            initState.(an).iEl = 0*onevec;
            initState.(an).i   = 0*onevec;

            nc = model.(elyte).G.getNumberOfCells();
            initState.(elyte).pi  = zeros(nc, 1);
            initState.(elyte).phi = zeros(nc, 1);

            onevec = ones(model.(ct).N, 1);

            initState.(ct).pi  = 0*onevec;
            initState.(ct).iEl = 0*onevec;
            initState.(ct).iHp = 0*onevec;
            initState.(ct).i   = 0*onevec;

            initState.(ctrl).I = 0;

            initState = model.evalVarName(initState, {an, 'Eocv'});
            initState.(ctrl).U = initState.(an).Eocv(1);

            initState.time = 0;

        end

        function model = setupForSimulation(model)

            model.isRootSimulationModel = true;

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
            forces.src   = [];
            
        end

        function model = validateModel(model, varargin)
        % do nothing
        end
        
    end

end
