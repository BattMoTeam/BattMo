classdef ProtonicMembraneGasSupply < BaseModel
    
    properties

        GasSupplyBc
        Control
        
        constants
        
        molecularWeights
        permeability
        viscosity
        T % Temperature
        control % Control structure
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        couplingTerms

        helpers
        
    end
    
    methods
        
        function model = ProtonicMembraneGasSupply(inputparams)
            
            model = model@BaseModel();

            fdnames = {'G'               , ...
                       'molecularWeights', ...
                       'permeability'    , ...
                       'viscosity'       , ...
                       'control'         , ...
                       'couplingTerms'   , ...
                       'T'};
            
            model = dispatchParams(model, inputparams, fdnames);

            model.GasSupplyBc = ProtonicMembraneGasSupplyBc(inputparams);
            model.Control     = BaseModel();
            
            model.constants = PhysicalConstants();
            
            model.gasInd.H2O = 1;
            model.gasInd.O2  = 2;
            model.nGas = 2;

            model = model.setupControl();

        end


        function model = setupForSimulation(model)
            
            model.isRootSimulationModel = true;

            shortNames = {'1', 'H2O';
                          '2', 'O2';
                          'massConses', 'massCons';
                          'GasSupplyBc', 'bc';
                          'controlEquations', 'ctrleqs';
                          'boundaryEquations', 'bceqs'};
            
            model = model.equipModelForComputation('shortNames', shortNames);
            
        end
        
        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            nGas = model.nGas;
            nCoup = numel(model.couplingTerms);
            
            varnames = {};

            % Components volume fractions
            varnames{end + 1} = VarName({}, 'volumefractions', nGas);
            % Total pressure
            varnames{end + 1} = 'pressure';
            % Gas densities 
            varnames{end + 1} = 'density';
            % Mass accumulation terms
            varnames{end + 1} = VarName({}, 'massAccums', nGas);
            % Mass Flux terms
            varnames{end + 1} = VarName({}, 'massFluxes', nGas);
            % Mass Source terms
            varnames{end + 1} = VarName({}, 'massSources', nGas);
            % Mass Conservation equations
            varnames{end + 1} = VarName({}, 'massConses', nGas);
            %% Control values
            % Rate (total mass) at the control
            % NOTE : Rate > 0 when injecting (opposity convention as massFlux)
            varnames{end + 1} = {'Control', 'rate'};
            % total pressure at the control
            varnames{end + 1} = {'Control', 'pressure'};
            % Pressure equation for the control pressure
            varnames{end + 1} = {'Control', 'pressureEq'};
            % Rate equation for the control equation
            varnames{end + 1} = {'Control', 'rateEq'};
            % control equation 
            varnames{end + 1} = {'Control', 'setupEq'};
            
            model = model.registerVarNames(varnames);

            fn = @ProtonicMembraneGasSupply.updateVolumeFraction;
            inputvarnames = {VarName({}, 'volumefractions', nGas, 1)};
            outputvarname = VarName({}, 'volumefractions', nGas, 2);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn = @(model, state) ProtonicMembraneGasSupply.updateDensity(model, state);
            fn = {fn, @(prop) PropFunction.literalFunctionCallSetupFn(prop)};
            inputvarnames = {VarName({}, 'volumefractions', nGas), 'pressure'};
            model = model.registerPropFunction({'density', fn, inputvarnames});

            fn = @ProtonicMembraneGasSupply.updateBCequations;
            inputvarnames = {VarName({'GasSupplyBc'}, 'massFluxes', nGas)     , ...
                             VarName({'GasSupplyBc'}, 'volumefractions', nGas), ...
                             {'GasSupplyBc', 'pressure'}                      , ...
                             {'GasSupplyBc', 'density'}                       , ...
                             VarName({}, 'volumefractions', nGas)             , ...
                             'pressure'                                       , ...
                             'density'};
            outputvarname = VarName({'GasSupplyBc'}, 'boundaryEquations', nGas);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn = @ProtonicMembraneGasSupply.updateControlPressureEquation;
            inputvarnames = {{'GasSupplyBc', 'pressure'}, {'Control', 'pressure'}};
            model = model.registerPropFunction({{'Control', 'pressureEq'}, fn, inputvarnames});

            fn = @ProtonicMembraneGasSupply.updateControlSetup;
            inputvarnames = {{'Control', 'pressure'}, ...
                             {'Control', 'rate'}};
            model = model.registerPropFunction({{'Control', 'setupEq'}, fn, inputvarnames});


            fn = @ProtonicMembraneGasSupply.updateBcVolumeFraction;
            inputvarnames = {'pressure', ...
                             {'GasSupplyBc', 'pressure'}, ...
                             VarName({}, 'volumefractions', nGas, 1)};
            outputvarname = VarName({'GasSupplyBc'}, 'volumefractions', nGas, 1);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            
            for igas = 1 : nGas

                fn = @ProtonicMembraneGasSupply.updateRateEq;
                inputvarnames = {VarName({'GasSupplyBc'}, 'massFluxes', nGas, igas), ...
                                 VarName({'Control'}, 'rate')};
                model = model.registerPropFunction({{'Control', 'rateEq'}, fn, inputvarnames});

                fn = @ProtonicMembraneGasSupply.updateMassFluxes;
                inputvarnames = {VarName({}, 'volumefractions', nGas, igas), ...
                                 'density', 'pressure'};
                outputvarname = VarName({}, 'massFluxes', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
                fn = @ProtonicMembraneGasSupply.updateMassAccums;
                fn = {fn, @(prop) PropFunction.accumFuncCallSetupFn(prop)};
                inputvarnames = {VarName({}, 'volumefractions', nGas, igas), ...
                                 'density'};
                outputvarname = VarName({}, 'massAccums', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                                
                fn = @ProtonicMembraneGasSupply.updateMassConses;
                inputvarnames = {VarName({}, 'massAccums', nGas, igas), ...
                                 VarName({}, 'massFluxes', nGas, igas), ...
                                 VarName({}, 'massSources', nGas, igas)};
                outputvarname = VarName({}, 'massConses', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @ProtonicMembraneGasSupply.updateMassSources;
                inputvarnames = {VarName({'GasSupplyBc'}, 'massFluxes', nGas, igas)};
                outputvarname = VarName({}, 'massSources', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
            end
            
        end

        function state = updateBcVolumeFraction(model, state)

            h = model.helpers;
            
            map1    = h.coupToBcMap;
            map2    = h.volumefractionMap;
            bccells = h.bccells;
            
            vfBc = map1*map2*(h.volumefractionValues);

            p    = state.pressure(bccells);
            vf   = state.volumefractions{1}(bccells);
            pBc  = state.GasSupplyBc.pressure;

            v = (p - pBc);

            vf(v < 0) = vfBc(v < 0);

            state.GasSupplyBc.volumefractions{1} = vf;
            
        end
        
        function state = updateControlSetup(model, state)

            helpers = model.helpers;

            eqs = {};

            map = helpers.rateMap;
            val = helpers.rateValues;

            eqs{end + 1} = 1./val.*(map*state.Control.rate) - 1; % equations are scaled to one

            map = helpers.pressureMap;
            val = helpers.pressureValues;

            eqs{end + 1} = 1./val.*(map*state.Control.pressure) - 1; % equations are scaled to one

            state.Control.setupEq = vertcat(eqs{:});

        end

        function state = updateMassSources(model, state)

            nGas = model.nGas;
            bccells = model.helpers.bccells;
            
            for igas = 1 : nGas
                src = 0.*state.pressure; % hacky initialization to get AD
                src(bccells) = state.GasSupplyBc.massFluxes{igas};
                srcs{igas} = - src;
            end

            state.massSources = srcs;
            
        end
        
        function model = setupControl(model)

            % Two crontol type indexed 'pressure-composition' : 1, 'flux-composition' : 2
            typetbl.type = [1; 2];
            typetbl = IndexArray(typetbl);

            % Two values are given for each control. They are indexed by comp
            % 'pressure-composition' : 1) pressure 2) H2O volume fraction
            % 'flux-composition' : 1) rate 2) H2O volume fraction
            %
            comptbl.comp = [1; 2];
            comptbl = IndexArray(comptbl);

            couplingTerms = model.couplingTerms;
            ctrl = model.control;

            ctrlvals = [];
            % The control values ctrlvals are stored in comptypecouptbl
            comptypecouptbl = IndexArray([], 'fdnames', {'comp', 'type', 'coup'});

            function icoupterm = findCouplingTerm(ictrl)
                name = ctrl(ictrl).name;
                icoupterm = [];
                for icoup = 1 : numel(couplingTerms)
                    if strcmp(name, couplingTerms{icoup}.name)
                        icoupterm = icoup;
                        return
                    end
                end
                error('couplingTerm structure not found');
            end
            
            for ictrl = 1 : numel(ctrl)

                ctrlinput = ctrl(ictrl);
                
                clear typetbl2
                switch ctrlinput.type
                  case 'pressure-composition'
                    typetbl2.type = 1;
                  case 'rate-composition'
                    typetbl2.type = 2;
                  otherwise
                    error('ctrlinput type not recognized');
                end
                typetbl2 = IndexArray(typetbl2);
                comptypecouptbl2 = crossIndexArray(comptbl, typetbl2, {});

                icoupterm = findCouplingTerm(ictrl);
                
                comptypecouptbl2 =  comptypecouptbl2.addInd('coup', icoupterm);
                vals2 = ctrlinput.values;

                comptypecouptbl = concatIndexArray(comptypecouptbl, comptypecouptbl2, {});
                ctrlvals = [ctrlvals; vals2];
                
            end

            couptbl.coup = (1 : numel(couplingTerms))';
            couptbl = IndexArray(couptbl);
            
            bccellfacecouptbl = IndexArray([], 'fdnames', {'faces', 'cells', 'coup'});

            for icoup = 1 : numel(couplingTerms)

                clear bccellfacecouptbl2;
                bccellfacecouptbl2.faces = couplingTerms{icoup}.couplingfaces;
                bccellfacecouptbl2.cells = couplingTerms{icoup}.couplingcells;
                bccellfacecouptbl2 = IndexArray(bccellfacecouptbl2);
                nc = bccellfacecouptbl2.num;
                bccellfacecouptbl2 =  bccellfacecouptbl2.addInd('coup', icoup*ones(nc, 1));

                bccellfacecouptbl = concatIndexArray(bccellfacecouptbl, bccellfacecouptbl2, {});

            end

            map          = TensorMap;
            map.fromTbl  = bccellfacecouptbl;
            map.toTbl    = couptbl;
            map.mergefds = {'coup'};

            M = SparseTensor;
            M = M.setFromTensorMap(map);
            M = M.getMatrix();

            helpers.bcToCoupMap = M;
            helpers.coupToBcMap = M';
            
            bccellfacecomptypecouptbl = crossIndexArray(bccellfacecouptbl, comptypecouptbl, {'coup'});

            %%  setup pressureMap
            
            clear comptypetbl2
            comptypetbl2.type = 1; % pressure index is 1
            comptypetbl2.comp = 1; % first index gives pressure value
            comptypetbl2 = IndexArray(comptypetbl2);
            comptypecouptbl2  = crossIndexArray(comptypetbl2, comptypecouptbl, {'type', 'comp'});
            
            map          = TensorMap();
            map.fromTbl  = couptbl;
            map.toTbl    = comptypecouptbl2;
            map.mergefds = {'coup'};

            M = SparseTensor();
            M = M.setFromTensorMap(map);
            M = M.getMatrix();

            helpers.pressureMap = M;

            map          = TensorMap();
            map.fromTbl  = comptypecouptbl;
            map.toTbl    = comptypecouptbl2;
            map.mergefds = {'coup', 'comp', 'type'};
            map          = map.setup();

            helpers.pressureValues = map.eval(ctrlvals);

            %%  setup rateMap
            
            clear comptypetbl2
            comptypetbl2.type = 2; % rate index is 2
            comptypetbl2.comp = 1; % first index gives rate value
            comptypetbl2 = IndexArray(comptypetbl2);
            comptypecouptbl2 = crossIndexArray(comptypetbl2, comptypecouptbl, {'type', 'comp'});
            
            map          = TensorMap();
            map.fromTbl  = couptbl;
            map.toTbl    = comptypecouptbl2;
            map.mergefds = {'coup'};

            M = SparseTensor();
            M = M.setFromTensorMap(map);
            M = M.getMatrix();

            helpers.rateMap = M;

            map          = TensorMap();
            map.fromTbl  = comptypecouptbl;
            map.toTbl    = comptypecouptbl2;
            map.mergefds = {'coup', 'comp', 'type'};
            map          = map.setup();

            helpers.rateValues = map.eval(ctrlvals);


            %%  setup volumeFractionMap
            
            clear comptbl2
            comptbl2.comp = 2; % first index gives rate value
            comptbl2 = IndexArray(comptbl2);
            comptypecouptbl2 = crossIndexArray(comptbl2, comptypecouptbl, {'comp'});
            
            map          = TensorMap();
            map.fromTbl  = couptbl;
            map.toTbl    = comptypecouptbl2;
            map.mergefds = {'coup'};

            M = SparseTensor();
            M = M.setFromTensorMap(map);
            M = M.getMatrix();

            helpers.volumefractionMap = M;

            map = TensorMap();
            map.fromTbl = comptypecouptbl;
            map.toTbl = comptypecouptbl2;
            map.mergefds = {'coup', 'comp', 'type'};
            map = map.setup();

            helpers.volumefractionValues = map.eval(ctrlvals);

            helpers.bccells           = bccellfacecouptbl.get('cells');
            helpers.bcfaces           = bccellfacecouptbl.get('faces');
            helpers.bccellfacecouptbl = bccellfacecouptbl;
            
            model.helpers = helpers;
            
        end


        function state = updateVolumeFraction(model, state)

            state.volumefractions{2} = 1 - state.volumefractions{1};
            
        end

        function state = updateControlPressureEquation(model, state)

            map = model.helpers.coupToBcMap;
            
            pcontrol = state.Control.pressure;
            pbc      = state.GasSupplyBc.pressure;

            state.Control.pressureEq = pbc - map*pcontrol;
            
        end


        function state = updateRateEq(model, state)

            map = model.helpers.bcToCoupMap;
            nGas = model.nGas;

            eq = state.Control.rate;
            
            for igas = 1 : nGas
                % Note sign due to convention
                eq = eq +  map*state.GasSupplyBc.massFluxes{igas};

            end

            state.Control.rateEq = eq;
            
        end

        function state = updateBCequations(model, state)

            nGas = model.nGas;
            K    = model.permeability;
            mu   = model.viscosity;
            
            bccells = model.helpers.bccells;
            bcfaces = model.helpers.bcfaces;

            Tbc = model.G.getBcTrans(bcfaces);

            pIn   = state.pressure(bccells);
            pBc   = state.GasSupplyBc.pressure;
            rhoBc = state.GasSupplyBc.density ; % note that this has been already "upwinded" (see specific update of volumefractions)
            vfsBc = state.GasSupplyBc.volumefractions; % note that those has been already "upwinded" (see specific update)

            for igas = 1 : nGas

                vfBc   = vfsBc{igas};
                bcFlux = state.GasSupplyBc.massFluxes{igas};

                bceqs{igas} = rhoBc.*vfBc*K/mu.*Tbc.*(pIn - pBc) - bcFlux;

            end

            state.GasSupplyBc.boundaryEquations = bceqs;
            
        end
        

        function state = updateMassFluxes(model, state)

            K    = model.permeability;
            mu   = model.viscosity;
            nGas = model.nGas;

            p   = state.pressure;
            vfs = state.volumefractions;
            rho = state.density;
            
            v = assembleFlux(model, p, rho.*K/mu);
            
            for igas = 1 : nGas

                state.massFluxes{igas} = assembleUpwindFlux(model, v, vfs{igas});
                
            end
            
        end
        
        function state = updateMassAccums(model, state, state0, dt)

            vols = model.G.getVolumes();
            nGas = model.nGas;

            vfs  = state.volumefractions;
            rho  = state.density;
            vfs0 = state0.volumefractions;
            rho0 = state0.density;
            
            for igas = 1 : nGas

                state.massAccums{igas} = vols.*(rho.*vfs{igas} - rho0.*vfs0{igas})/dt;
                
            end
            
        end

        
        function state = updateMassConses(model, state)

            nGas = model.nGas;

            for igas = 1 : nGas

                state.massConses{igas} = state.massAccums{igas} + model.G.getDiv(state.massFluxes{igas}) - state.massSources{igas};
                
            end
            
        end
        
        function initstate = setupInitialState(model)

            nGas    = model.nGas;
            nc      = model.G.getNumberOfCells();
            nbc     = numel(model.helpers.bcfaces);
            gasInd  = model.gasInd;
            control = model.control;


            names = arrayfun(@(ctrl) ctrl.name, control, 'uniformoutput', false);
            [lia, locb] = ismember('External output', names);
            assert(lia, 'External output control not found');
            control = control(locb);

            assert(strcmp(control.type, 'pressure-composition'), 'Here we expect pressure controled output');
            
            values = control.values;

            p  = values(1);
            vf = values(2);
            
            initstate.pressure           = p*ones(nc, 1);
            initstate.volumefractions{1} = vf *ones(nc, 1);
            initstate.volumefractions{2} = (1 - vf) *ones(nc, 1);
            
            initstate.GasSupplyBc.pressure               = p*ones(nbc, 1);
            initstate.GasSupplyBc.massFluxes{gasInd.H2O} = zeros(nbc, 1);
            initstate.GasSupplyBc.massFluxes{gasInd.O2}  = zeros(nbc, 1);
            initstate.GasSupplyBc.volumefractions{1}     = vf*ones(nbc, 1);
            
            nctrl = numel(model.control);
            initstate.Control.rate               = zeros(nctrl, 1);
            initstate.Control.pressure           = p*ones(nctrl, 1);
            initstate.Control.volumefractions{1} = vf*ones(nctrl, 1);

            initstate = model.evalVarName(initstate, VarName({}, 'volumefractions', nGas, 2));
            initstate = model.evalVarName(initstate, VarName({}, 'density'));
            
        end

        function cleanState = addStaticVariables(model, cleanState, state, state0)
            
            cleanState = addStaticVariables@BaseModel(model, cleanState, state);

            cleanState.volumefractions{2} = state.volumefractions{2};
            cleanState.density            = state.density;
            
        end

        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@PhysicalModel(model);
            forces.src = [];
            
        end

        function model = validateModel(model, varargin)
        % do nothing
        end
        
    end

    methods(Static)

        function state = updateDensity(model, state)

            Mws  = model.molecularWeights;
            T    = model.T;
            nGas = model.nGas;
            c    = model.constants;

            p   = state.pressure;
            vfs = state.volumefractions;

            density = 0*p;
            
            for igas = 1 : nGas

                density = density + (p.*vfs{igas}.*Mws(igas))./(c.R*T);
                
            end

            state.density = density;
            
        end

    end
    
end
