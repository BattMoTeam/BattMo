classdef ProtonicMembraneGasSupply < BaseModel
    
    properties

        GasSupplyBc
        
        constants
        
        molecularWeights
        permeability
        viscosity
        T % Temperature
        control % Control structure
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        couplingTerms
        
    end
    
    methods
        
        function model = ProtonicMembraneGasSupply(paramobj)
            
            model = model@BaseModel();

            fdnames = {'G'               , ...
                       'molecularWeights', ...
                       'permeability'    , ...
                       'viscosity'       , ...
                       'control'         , ...
                       'couplingTerms'   , ...
                       'T'};
            
            model = dispatchParams(model, paramobj, fdnames);

            model.GasSupplyBc = ProtonicMembraneGasSupplyBc([]);
            
            model.GasSupplyBc.molecularWeights = model.molecularWeights;
            model.GasSupplyBc.T                = model.T;
            
            model.operators = localSetupOperators(model.G);
            
            model.constants = PhysicalConstants();
            
            model.gasInd.H2O = 1;
            model.gasInd.O2  = 2;
            model.nGas = 2;

            model = model.setupControl();

        end


        function model = setupForSimulation(model)
            
            model.isSimulationModel = true;

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


            for igas = 1 : nGas

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

        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@PhysicalModel(model);
            forces.src = [];
            
        end

        function state = updateMassSources(model, state)

            nGas = model.nGas;
            bccells = model.GasSupplyBc.controlHelpers.bccells;
            
            for igas = 1 : nGas
                src = 0.*state.pressure; % hacky initialization to get AD
                src(bccells) = state.GasSupplyBc.massFluxes{igas};
                srcs{igas} = - src;
            end

            state.massSources = srcs;
            
        end
        
        function model = setupControl(model)

            % Two components indexed H2O : 1, O2 : 2
            comptbl.comp = [1; 2];
            comptbl = IndexArray(comptbl);
            
            % Two crontol type indexed 'pressure' : 1, 'flux' : 2
            typetbl.type = [1; 2];
            typetbl = IndexArray(typetbl);

            couplingTerms = model.couplingTerms;
            ctrl = model.control;

            ctrlvals = [];
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
                  case 'flux'
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

            bccellfacecomptypecouptbl = crossIndexArray(bccellfacecouptbl, comptypecouptbl, {'coup'});

            for itype = 1 : typetbl.num
                
                for icomp = 1 : comptbl.num
                    
                    clear comptypetbl
                    comptypetbl.type = itype;
                    comptypetbl.comp = icomp;
                    comptypetbl = IndexArray(comptypetbl);
                    
                    comptypecouptbl2 = crossIndexArray(comptypecouptbl, comptypetbl, {'comp', 'type'});

                    bccellfacecomptypecouptbl = crossIndexArray(comptypecouptbl2, bccellfacecouptbl, {'coup'});
                
                    map = TensorMap();
                    map.fromTbl = bccellfacecouptbl;
                    map.toTbl = bccellfacecomptypecouptbl;
                    map.mergefds = {'cells', 'faces', 'coup'};
                    map = map.setup();

                    controlMap = SparseTensor();
                    controlMap = controlMap.setFromTensorMap(map);
                    controlMap = controlMap.getMatrix();
                    
                    controlMaps{itype, icomp} = controlMap;

                    map = TensorMap();
                    map.fromTbl = comptypecouptbl;
                    map.toTbl = bccellfacecomptypecouptbl;
                    map.mergefds = {'comp', 'type', 'coup'};
                    map = map.setup();

                    controlVals{itype, icomp} = map.eval(ctrlvals);
                    
                end
                
            end

            controlHelpers.maps              = controlMaps;
            controlHelpers.vals              = controlVals;
            controlHelpers.bccells           = bccellfacecouptbl.get('cells');
            controlHelpers.bcfaces           = bccellfacecouptbl.get('faces');
            controlHelpers.bccellfacecouptbl = bccellfacecouptbl;
            
            model.GasSupplyBc.controlHelpers = controlHelpers;
            
        end


        function state = updateVolumeFraction(model, state)

            state.volumefractions{2} = 1 - state.volumefractions{1};
            
        end
        
        function initstate = setupInitialState(model)

            nc      = model.G.cells.num;
            nbc     = model.GasSupplyBc.getNumberBcFaces();
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
            
            initstate.GasSupplyBc.pressure               = p*ones(nbc, 1);
            initstate.GasSupplyBc.volumefractions{1}     = vf *ones(nbc, 1);
            initstate.GasSupplyBc.massFluxes{gasInd.H2O} = zeros(nbc, 1);
            initstate.GasSupplyBc.massFluxes{gasInd.O2}  = zeros(nbc, 1);
            
        end

        
        function state = updateBCequations(model, state)

            nGas = model.nGas;
            K    = model.permeability;
            mu   = model.viscosity;
            
            bccells = model.GasSupplyBc.controlHelpers.bccells;
            bcfaces = model.GasSupplyBc.controlHelpers.bcfaces;

            Tbc = model.operators.transFaceBC(bcfaces);
            Tbc = K/mu.*Tbc;

            pIn   = state.pressure(bccells);
            pBc   = state.GasSupplyBc.pressure;
            rhoIn = state.density(bccells);
            rhoBc = state.GasSupplyBc.density;
            vfsIn = state.volumefractions;
            vfsBc = state.GasSupplyBc.volumefractions;

            for igas = 1 : nGas

                vfIn   = vfsIn{igas}(bccells);
                vfBc   = vfsBc{igas};
                bcFlux = state.GasSupplyBc.massFluxes{igas};
                
                bceqs{igas} = ProtonicMembraneGasSupplyBc.setupBcEquation(model , ...
                                                                          bcFlux, ...
                                                                          pIn   , ...
                                                                          pBc   , ...
                                                                          rhoIn , ...
                                                                          rhoBc , ...
                                                                          vfIn  , ...
                                                                          vfBc  , ...
                                                                          Tbc);

            end

            state.GasSupplyBc.boundaryEquations = bceqs;
            
        end
        
        function cleanState = addStaticVariables(model, cleanState, state, state0)
            
            cleanState = addStaticVariables@BaseModel(model, cleanState, state);

            cleanState.volumefractions{2} = state.volumefractions{2};
            cleanState.density            = state.density;
            
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

            vols = model.G.cells.volumes;
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
            op   = model.operators;

            for igas = 1 : nGas

                state.massConses{igas} = state.massAccums{igas} + op.Div(state.massFluxes{igas}) - state.massSources{igas};
                
            end
            
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
