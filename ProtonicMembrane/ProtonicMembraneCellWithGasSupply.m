classdef ProtonicMembraneCellWithGasSupply < BaseModel
    
    properties
        
        % Temperature
        T
        % Structure with physical constants
        constants

        Cell
        GasSupply
        Interface
        
        couplingTerm

        helpers
        pmin  = 1e-3*barsa % cutoff value used in updateState method
        vfmin = 1e-8       % cutoff value used in updateState method
        vfmax = 1 - 1e-8   % cutoff value used in updateState method
        
    end
    
    methods
        
        function model = ProtonicMembraneCellWithGasSupply(paramobj)

            model = model@BaseModel();

            fdnames = {'T', ...
                       'couplingTerm'};
            
            model = dispatchParams(model, paramobj, fdnames);

            model.Cell      = ProtonicMembraneCell(paramobj.Cell);
            model.GasSupply = ProtonicMembraneGasSupply(paramobj.GasSupply);
            model.Interface = BaseModel();
            
            model = model.setupInterface();
            
            model.constants = PhysicalConstants();
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            gs    = 'GasSupply';
            ce    = 'Cell';
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';

            nGas = model.GasSupply.nGas;

            varnames = {};

            % Coupling equations
            varnames{end + 1} = VarName({}, 'massCouplingEquations', nGas);
            % We add time variables to this model and the Cell model. Needed for the control
            varnames{end + 1} = 'time';
            % Coupling damping variable
            varnames{end + 1} = 'beta';
            % pressure variable at the interface
            varnames{end + 1} = VarName({'Interface'}, 'pressures', nGas);
            % massFluxes variable at the interface
            varnames{end + 1} = VarName({'Interface'}, 'massFluxes', nGas);
            % massFlux equation at Interface (used to decouple)
            varnames{end + 1} = VarName({'Interface'}, 'massFluxEqs', nGas);

            
            model = model.registerVarNames(varnames);

            fn =  @ProtonicMembraneCellWithGasSupply.updateInterfaceEquation;
            inputvarnames = {{'GasSupply', 'pressure'}                      , ...
                             VarName({'GasSupply'}, 'volumefractions', nGas), ...
                             VarName({'Interface'}, 'pressures', nGas)      , ...
                             VarName({'Interface'}, 'massFluxes', nGas)};
            outputvarname = VarName({'Interface'}, 'massFluxEqs', nGas);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn =  @ProtonicMembraneCellWithGasSupply.updateAnodePressures;
            inputvarnames = {VarName({'Interface'}, 'pressures', nGas)};
            outputvarname = VarName({'Cell', 'Anode'}, 'pressures', nGas);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn =  @ProtonicMembraneCellWithGasSupply.updateCouplingEquation;
            inputvarnames = {{'Cell', 'Anode' 'iHp'}, ...
                             VarName({'Interface'}, 'massFluxes', nGas), ...
                             'beta'};
            outputvarname = VarName({}, 'massCouplingEquations', nGas);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            
            fn =  @ProtonicMembraneCellWithGasSupply.updateGasSupplySource;
            inputvarnames = {VarName({'Interface'}, 'massFluxes', nGas), VarName({'GasSupply', 'GasSupplyBc'}, 'massFluxes', nGas)};
            outputvarname = VarName({'GasSupply'}, 'massSources', nGas);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            
            model = model.setAsStaticVarName('time');
            
            fn =  @ProtonicMembraneCellWithGasSupply.dispatchTime;
            model = model.registerPropFunction({{ce, 'time'}, fn, {'time'}});

            fn = @ProtonicMembraneCell.updateBeta;
            fn = {fn, @(propfunction) PropFunction.drivingForceFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({'beta', fn, {'time'}});
            
        end

        function initstate = setupInitialState(model)

            gs = 'GasSupply';
            ce = 'Cell';
            an = 'Anode';
            
            gasInd = model.(gs).gasInd;
            nGas   = model.(gs).nGas;
            
            initstate.(gs) = model.(gs).setupInitialState();

            p = initstate.(gs).pressure;
            s = initstate.(gs).volumefractions{1};

            initstate.(gs).volumefractions{2} = 1 - initstate.(gs).volumefractions{1};

            for igas = 1 : nGas
                pgs{igas} = initstate.(gs).volumefractions{igas}.*initstate.(gs).pressure;
            end

            pH2O = pgs{gasInd.H2O};
            pO2  = pgs{gasInd.O2};

            model.(ce).(an).pH2O = pH2O(1);
            model.(ce).(an).pO2  = pO2(1);
            
            model.(ce) = model.(ce).setupComputationalGraph();

            initstate.(ce) = model.(ce).setupInitialState();

            initstate = model.evalVarName(initstate, VarName({gs}, 'density'));
            initstate.time = 0;

            map = model.helpers.interfaceMap;
            
            for igas = 1 : nGas

                initstate.Interface.pressures{igas} = map*pgs{igas};
                initstate.Interface.massFluxes{igas} = map*(0*pgs{igas});
                
            end
            
        end

        function model = setupInterface(model)

            coupterms = model.Cell.couplingTerms;
            
            for icoup = 1 : numel(coupterms)
                coupterm = coupterms{icoup};
                [lia, locb] = ismember('Anode', coupterm.componentnames);
                if lia
                    loc = true(2, 1);
                    loc(locb) = false;
                    faces = coupterm.couplingfaces(:, loc);
                    break
                end
            end

            faces = model.Cell.Electrolyte.G.mappings.parentGrid.mappings.facemap(faces);

            G = model.GasSupply.G;

            invfacemap = zeros(G.mappings.parentGrid.faces.num, 1);
            invfacemap(G.mappings.facemap) = (1 : G.faces.num)';

            faces = invfacemap(faces);

            tbls = setupSimpleTables(G);

            cellfacetbl = tbls.cellfacetbl;
            celltbl     = tbls.celltbl;
            
            bcfacetbl.faces = faces;
            bcfacetbl = IndexArray(bcfacetbl);
            bcfaceindtbl = bcfacetbl.addInd('ind', (1 : bcfacetbl.num)');

            bccellfaceindtbl = crossIndexArray(bcfaceindtbl, cellfacetbl, {'faces'});
            bccellfaceindtbl = sortIndexArray(bccellfaceindtbl, {'ind', 'faces', 'cells'});

            map = TensorMap();
            map.fromTbl = celltbl;
            map.toTbl = bccellfaceindtbl;
            map.mergefds = {'cells'};

            M = SparseTensor();
            M = M.setFromTensorMap(map);
            M = M.getMatrix();

            T = model.GasSupply.operators.T_all;
            T = T(bccellfaceindtbl.get('faces'));
            
            model.helpers.interfaceMap   = M;
            model.helpers.interfaceT     = T;
            model.helpers.interfaceCells = bccellfaceindtbl.get('cells');
            model.helpers.interfaceFaces = bccellfaceindtbl.get('faces');
            
        end

        function state = updateBeta(model, state, drivingForces)
            
            time = state.time;
            [~, ~, beta] = drivingForces.src(time);
            state.beta = beta;
            
        end

        function state = dispatchTime(model, state)

            state.Cell.time = state.time;
            
        end

        function state = updateInterfaceEquation(model, state)

            map  = model.helpers.interfaceMap;
            T    = model.helpers.interfaceT;
            nGas = model.GasSupply.nGas;
            K    = model.GasSupply.permeability;
            mu   = model.GasSupply.viscosity;
            
            pIn     = state.GasSupply.pressure;
            vfIns   = state.GasSupply.volumefractions;
            rho     = state.GasSupply.density;
            pBcs    = state.Interface.pressures;
            mfluxes = state.Interface.massFluxes;

            rho = map*rho;
            pIn = map*pIn;
            pBc = 0*pBcs{1}; % hacky initialisation
            mflux = 0*pBcs{1}; % hacky initialisation
            for igas = 1 : nGas
                vfIns{igas} = map*vfIns{igas};
                pBc = pBc + pBcs{igas};
                mflux = mflux + mfluxes{igas};
            end

            eqs{1} = mflux - rho.*K/mu.*T.*(pIn - pBc);
            eqs{2} = vfIns{1}.*pBc - pBcs{1};

            state.Interface.massFluxEqs = eqs;
            
        end


        function state = updateGasSupplySource(model, state)

            nGas     = model.GasSupply.nGas;
            bccells  = model.GasSupply.GasSupplyBc.controlHelpers.bccells;
            intcells = model.helpers.interfaceCells;
            
            for igas = 1 : nGas

                src = 0.*state.GasSupply.pressure; % hacky initialization to get AD

                src(bccells)  = src(bccells) + state.GasSupply.GasSupplyBc.massFluxes{igas};
                src(intcells) = src(intcells) + state.Interface.massFluxes{igas};

                srcs{igas} = - src;
                
            end

            state.GasSupply.massSources = srcs;
            
        end

        
        function state = updateCouplingEquation(model, state)

            nGas   = model.GasSupply.GasSupplyBc.nGas;
            gasInd = model.GasSupply.GasSupplyBc.gasInd;
            mws    = model.GasSupply.molecularWeights;
            F      = PhysicalConstants.F;
            
            iHp     = state.Cell.Anode.iHp;
            mfluxes = state.Interface.massFluxes;
            beta    = state.beta;
            
            % Chemical equation : 1/2*H2O <-> H^+ + e^- + 1/4*O2
            % Flux in the boundary conditions are oritented outwards from GasSupply domain.
            igas = gasInd.H2O;
            coupeqs{igas} = mfluxes{igas} - beta*1/2*mws(igas)*iHp/F;

            igas = gasInd.O2;
            coupeqs{igas} = mfluxes{igas} + beta*1/4*mws(igas)*iHp/F;

            state.massCouplingEquations = coupeqs;
            
        end

        function state = updateAnodePressures(model, state)

            state.Cell.Anode.pressures = state.Interface.pressures;
            
        end


        function model = setupForSimulation(model)

            model.isSimulationModel = true;

            shortNames = {'1'                    , 'H2O';
                          '2'                    , 'O2';
                          'massConses'           , 'massCons';
                          'Cell'                 , 'ce';
                          'GasSupply'            , 'gs';
                          'GasSupplyBc'          , 'gsBc';
                          'massCouplingEquations', 'massCoupEqs';
                          'controlEquations'     , 'ctrlEqs';
                          'boundaryEquations'    , 'bcEqs';
                          'Electrolyte'          , 'elyte';
                          'Anode'                , 'an';
                          'Cathode'              , 'ct';
                          'Control'              , 'ctrl'};

            model = model.equipModelForComputation('shortNames', shortNames);

        end

        
        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@BaseModel(model);
            forces.src   = [];
            
        end

        function model = validateModel(model, varargin)
        % do nothing
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            
            [state, report] = updateState@BaseModel(model, state, problem, dx, drivingForces);

            % cutoff the pressures

            gs    = 'GasSupply';
            ce    = 'Cell';
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';

            nGas = model.(gs).nGas;

            pmin  = model.pmin;
            vfmin = model.vfmin;
            vfmax = model.vfmax;
            
            state.(gs).pressure             = max(pmin, state.(gs).pressure);
            state.(gs).GasSupplyBc.pressure = max(pmin, state.(gs).GasSupplyBc.pressure);
            
            state.(gs).volumefractions{1}             = max(vfmin, state.(gs).volumefractions{1});
            state.(gs).GasSupplyBc.volumefractions{1} = max(vfmin, state.(gs).GasSupplyBc.volumefractions{1});
            
            state.(gs).volumefractions{1}             = min(vfmax, state.(gs).volumefractions{1});
            state.(gs).GasSupplyBc.volumefractions{1} = min(vfmax, state.(gs).GasSupplyBc.volumefractions{1});
            
        end
        
    end
    
end
