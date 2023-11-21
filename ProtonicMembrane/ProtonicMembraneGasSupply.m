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
        
        standalone
        
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

            model.operators = localSetupOperators(model.G);
            
            model.constants = PhysicalConstants();
            
            model.gasInd.H2O = 1;
            model.gasInd.O2  = 2;
            model.nGas = 2;

            model = model.setupControl();
            
        end

        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);

            nGas = model.nGas;
            
            varnames = {};

            % Gas partial pressures
            varnames{end + 1} = VarName({}, 'pressures', nGas);
            % Total pressure
            varnames{end + 1} = 'pressure';
            % Gas densities 
            varnames{end + 1} = VarName({}, 'densities', nGas);
            % Mass accumulation terms
            varnames{end + 1} = VarName({}, 'massAccums', nGas);
            % Mass Flux terms
            varnames{end + 1} = VarName({}, 'massFluxes', nGas);
            % Mass Source terms
            varnames{end + 1} = VarName({}, 'massSources', nGas);
            % Mass Conservation equations
            varnames{end + 1} = VarName({}, 'massConses', nGas);

            model = model.registerVarNames(varnames);

            fn = @ProtonicMembraneGasSupply.updatePressure;
            inputvarnames = {VarName({}, 'pressures', nGas)};
            model = model.registerPropFunction({'pressure', fn, inputvarnames})
            
            for igas = 1 : nGas
                
                fn = @(model, state) ProtonicMembraneGasSupply.updateDensities(model, state);
                fn = {fn, @(prop) PropFunction.literalFunctionCallSetupFn(prop)};
                inputvarnames = {VarName({}, 'pressures', nGas, igas)};
                outputvarname = VarName({}, 'densities', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @ProtonicMembraneGasSupply.massFluxes;
                inputvarnames = {VarName({}, 'pressures', nGas, igas), 'pressure'};
                outputvarname = VarName({}, 'massFluxes', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
                fn = @ProtonicMembraneGasSupply.updateMassAccums;
                fn = {fn, @(prop) PropFunction.accumFuncCallSetupFn(prop)};
                inputvarnames = {VarName({}, 'densities', nGas, igas)};
                outputvarname = VarName({}, 'massAccums', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                                
                fn = @ProtonicMembraneGasSupply.updateMassConses;
                inputvarnames = {VarName({}, 'massAccums', nGas, igas)  , ...
                                 VarName({}, 'massFluxes', nGas, igas) , ...
                                 VarName({}, 'massSources', nGas, igas)};
                outputvarname = VarName({}, 'massConses', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @ProtonicMembraneGasSupply.updateBcFluxTerms;
                inputvarnames = {VarName({'GasSupplyBc'}, 'massFluxes', nGas, igas), ...
                                 VarName({'GasSupplyBc'}, 'densities', nGas, igas) , ...
                                 VarName({'GasSupplyBc'}, 'pressures', nGas, igas) , ...
                                 VarName({}, 'densities', nGas, igas)              , ...
                                 VarName({}, 'pressures', nGas, igas)};
                outputvarname = VarName({'GasSupplyBc'}, 'boundaryEquations', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @ProtonicMembraneGasSupply.updateMassSources;
                inputvarnames = {VarName({'GasSupplyBc'}, 'massFluxes', nGas, igas)};
                outputvarname = VarName({}, 'massSources', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                
            end
        end

        function model = setupControl(model)

            % Two components indexed H2O : 1, O2 : 2
            comptbl.comp = [1; 2];
            comptbl = IndexArray(comptbl);
            
            % Two crontol type indexed 'pressure' : 1, 'flux' : 2
            typetbl.type = [1; 2];
            typetbl = IndexArray(typetbl);

            couplingTerms = model.couplingTerms;
            ctrl          = model.control;

            assert(numel(couplingTerms) == numel(ctrl), 'mismatch between coupling terms input');

            ctrlvals = [];
            comptypecouptbl = IndexArray([], 'fdnames', {'comp', 'type', 'coup'});

            for icoup = 1 : numel(ctrl)

                ctrlinput = ctrl(icoup);
                
                clear typetbl2
                switch ctrlinput.type
                  case 'pressure'
                    typetbl2.type = 1;
                  case 'flux'
                    typetbl2.type = 2;
                  otherwise
                    error('ctrlinput type not recognized');
                end
                typetbl2 = IndexArray(typetbl2);
                comptypecouptbl2 = crossIndexArray(comptbl, typetbl2, {});
                comptypecouptbl2 =  comptypecouptbl2.addInd('coup', icoup);
                vals2 = ctrlinput.values;

                comptypecouptbl = concatIndexArray(comptypecouptbl, comptypecouptbl2, {});
                ctrlvals = [ctrlvals; vals2];
            end

            bccellfacecouptbl.faces = [];
            bccellfacecouptbl.cells = [];
            bccellfacecouptbl.coup = [];
            bccellfacecouptbl = IndexArray(bccellfacecouptbl);
            
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

            bccellfacetbl = projIndexArray(bccellfacecomptypecouptbl, {'cells', 'faces'});

            for itype = 1 : typetbl.num
                
                for icomp = 1 : comptbl.num
                    
                    clear comptypetbl
                    comptypetbl.type = itype;
                    comptypetbl.comp = icomp;
                    comptypetbl = IndexArray(comptypetbl);
                    
                    comptypecouptbl2 = crossIndexArray(comptypecouptbl, comptypetbl, {'comp', 'type'});

                    bccellfacecomptypecouptbl = crossIndexArray(comptypecouptbl2, bccellfacecouptbl, {'coup'});
                
                    map = TensorMap();
                    map.fromTbl = bccellfacetbl;
                    map.toTbl = bccellfacecomptypecouptbl;
                    map.mergefds = {'cells', 'faces'};
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

            controlHelpers.maps    = controlMaps;
            controlHelpers.vals    = controlVals;
            controlHelpers.bccells = bccellfacetbl.get('cells');
            controlHelpers.bcfaces = bccellfacetbl.get('faces');
            
            model.GasSupplyBc.controlHelpers = controlHelpers;
            
        end
        
        function state = updatePressure(model, state)

            nGas = model.nGas;
            
            p = state.pressures{1};

            for igas = 2 : nGas
                p = p + state.pressures{2};
            end

            state.pressure = p;
            
        end
        

        function state = updateBcFluxTerms(model, state)

            nGas = model.nGas;
            K    = model.permeability;
            mu   = model.viscosity;
            
            bccells = model.GasSupplyBc.controlHelpers.bccells;
            bcfaces = model.GasSupplyBc.controlHelpers.bcfaces;

            Tbc = model.operators.transFaceBC(bcfaces);
            Tbc = K/mu.*Tbc;
            
            for igas = 1 : nGas

                pIn   = state.pressures{igas}(bccells);
                pBc   = state.GasSupplyBc.pressures{igas};
                rhoIn = state.densities{igas}(bccells);
                rhoBc = state.GasSupplyBc.densities{igas};

                bcFlux = state.GasSupplyBc.massFluxes{igas};
                
                bceqs{igas} = ProtonicMembraneGasSupplyBc.setupBcEquation(model , ...
                                                                         bcFlux, ...
                                                                         pIn   , ...
                                                                         pBc   , ...
                                                                         rhoIn , ...
                                                                         rhoBc , ...
                                                                         Tbc);

            end

            state.GasSupplyBc.boundaryEquations = bceqs;
            
        end
        
        
        function state = massFluxes(model, state)

            K    = model.permeability;
            mu   = model.viscosity;
            nGas = model.nGas;
            
            v = assembleHomogeneousFlux(model, state.pressure, K/mu);
            
            for igas = 1 : nGas

                rho = state.densities{igas};
                state.massFluxes{igas} = assembleUpwindFlux(model, v, rho);
                
            end
            
        end
        
        function state = updateMassAccums(model, state, state0, dt)

            vols = model.G.cells.volumes;
            nGas = model.nGas;
            
            for igas = 1 : nGas

                rho  = state.densities{igas};
                rho0 = state0.densities{igas};

                state.massAccums{igas} = vols.*(rho - rho0)/dt;
                
            end
            
        end

        
        function state = updateMassConses(model, state)

            nGas = model.nGas;
            op   = model.operators;

            for igas = 1 : nGas

                state.massConses{igas} = state.massAccums{igas} + op.div(state.massFluxes{igas}) + state.massSources{igas};
                
            end
            
        end

        
        
    end

    methods(Static)

        function state = updateDensities(model, state)

            Mws  = model.molecularWeights;
            T    = model.T;
            nGas = model.nGas;
            
            for igas = 1 : nGas

                state.densities{igas} = Mws(igas)*state.pressures{igas}.*c.R.*T;

            end
            
        end

    end
    
end
