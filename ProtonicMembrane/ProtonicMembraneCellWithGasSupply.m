classdef ProtonicMembraneCellWithGasSupply < BaseModel
    
    properties
        
        % Temperature
        T
        % Structure with physical constants
        constants

        Cell
        GasSupply

        couplingTerm

        helpers
        
    end
    
    methods
        
        function model = ProtonicMembraneCellWithGasSupply(paramobj)

            model = model@BaseModel();

            fdnames = {'T', ...
                       'couplingTerm'};
            
            model = dispatchParams(model, paramobj, fdnames);

            model.Cell = ProtonicMembraneCell(paramobj.Cell);

            % We setup coupling terms for the Gas Supply
            coupTerm = model.couplingTerm;
            couplingnames = cellfun(@(name) coupTerm.getComponentName(name), coupTerm.componentnames, 'uniformoutput', false);

            [lia, lib] = ismember('GasSupply', couplingnames);

            assert(lia, 'GasSupply coupling not found');
            
            coupTerm.name           = 'External Anode';
            coupTerm.componentnames = {'External'};
            coupTerm.couplingcells  = coupTerm.couplingcells(:, lib);
            coupTerm.couplingfaces  = coupTerm.couplingfaces(:, lib);

            paramobj.GasSupply.couplingTerms{end + 1} = coupTerm;
            ncoup = numel(paramobj.GasSupply.couplingTerms);

            model.GasSupply = ProtonicMembraneGasSupply(paramobj.GasSupply);
            
            model = setupGasSupplyCellCoupling(model, ncoup);
            
            model.constants = PhysicalConstants();
            
        end

        function initstate = setupInitialState(model)

            gs = 'GasSupply';
            ce = 'Cell';
            an = 'Anode';
            
            gasInd  = model.(gs).gasInd;

            initstate.(gs) = model.(gs).setupInitialState();
            
            pH2O = initstate.(gs).pressures{gasInd.H2O};
            pO2  = initstate.(gs).pressures{gasInd.O2};

            model.(ce).(an).pH2O = pH2O(1);
            model.(ce).(an).pO2  = pO2(1);
            
            model.(ce) = model.(ce).setupComputationalGraph();

            initstate.(ce) = model.(ce).setupInitialState();

            initstate = model.evalVarName(initstate, VarName({gs}, 'densities', 2));
            initstate.time = 0;

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
            varnames{end + 1} = {ce, 'time'};
            
            model = model.registerVarNames(varnames);

            fn =  @ProtonicMembraneCellWithGasSupply.updateBcAnodePressures;
            inputvarnames = {VarName({'GasSupply', 'GasSupplyBc'}, 'pressures', nGas)};
            outputvarname = VarName({'Cell', 'Anode'}, 'pressures', nGas);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn =  @ProtonicMembraneCellWithGasSupply.updateCouplingEquation;
            inputvarnames = {{'Cell', 'Anode' 'jHp'}, VarName({'GasSupply', 'GasSupplyBc'}, 'massFluxes', nGas)};
            outputvarname = VarName({}, 'massCouplingEquations', nGas);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            model = model.setAsStaticVarName('time');
            
            fn =  @ProtonicMembraneCellWithGasSupply.dispatchTime;
            model = model.registerPropFunction({{ce, 'time'}, fn, {'time'}});
            
        end

        function state = dispatchTime(model, state)

            state.Cell.time = state.time;
            
        end
            
        
        function state = updateCouplingEquation(model, state)

            nGas   = model.GasSupply.GasSupplyBc.nGas;
            gasInd = model.GasSupply.GasSupplyBc.gasInd;
            mws    = model.GasSupply.molecularWeights;
            map    = model.helpers.couplingMap;
            F      = PhysicalConstants.F;
            
            jHp      = state.Cell.Anode.jHp;
            mfluxbcs = state.GasSupply.GasSupplyBc.massFluxes;

            % Chemical equation : 1/2*H2O <-> H^+ + e^- + 1/4*O2
            % Flux in the boundary conditions are oritented outwards from GasSupply domain.
            igas = gasInd.H2O;
            coupeqs{igas} = map*mfluxbcs{igas} - 1/2*mws(igas)*F*jHp;

            igas = gasInd.O2;
            coupeqs{igas} = map*mfluxbcs{igas} + 1/4*mws(igas)*F*jHp;

            state.massCouplingEquations = coupeqs;
            
        end

        function state = updateBcAnodePressures(model, state)

            nGas = model.GasSupply.GasSupplyBc.nGas;
            map  = model.helpers.couplingMap;

            for igas = 1 : nGas

                state.Cell.Anode.pressures{igas} = map*state.GasSupply.GasSupplyBc.pressures{igas};
                
            end
            
        end

        
        function model = setupGasSupplyCellCoupling(model, ncoup)


            bccellfacecouptbl = model.GasSupply.GasSupplyBc.controlHelpers.bccellfacecouptbl;

            coupTerm = model.couplingTerm;

            componentnames = cellfun(@(name) coupTerm.getComponentName(name), coupTerm.componentnames, 'uniformoutput', false);

            [lia, locb] = ismember('GasSupply', componentnames);
            assert(lia, 'could not find GasSupply');

            anodecellfacecouptbl.cells = coupTerm.couplingcells(:, locb);
            anodecellfacecouptbl.faces = coupTerm.couplingfaces(:, locb);
            anodecellfacecouptbl = IndexArray(anodecellfacecouptbl);

            anodecellfacecouptbl = anodecellfacecouptbl.addInd('coup', ncoup);
            
            map = TensorMap();
            map.fromTbl = bccellfacecouptbl;
            map.toTbl = anodecellfacecouptbl;
            map.mergefds = {'cells', 'faces', 'coup'};

            M = SparseTensor();
            M = M.setFromTensorMap(map);
            M = M.getMatrix();

            model.helpers.couplingMap = M;
            
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
                          'controlEquation'      , 'ctrlEq';
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
        
        
    end
    
end
