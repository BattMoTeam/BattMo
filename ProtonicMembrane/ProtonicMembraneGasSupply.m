classdef ProtonicMembraneGasSupply < BaseModel
    
    properties

        GasSupplyBc
        
        constants
        
        molecularWeights
        permeability
        viscosity
        T % Temperature
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        coupTerms
        
        standalone
        
    end
    
    methods
        
        function model = ProtonicMembraneGasSupply(paramobj)
            
            model = model@BaseModel();

            fdnames = {'G'               , ...
                       'molecularWeights', ...
                       'permeability'    , ...
                       'viscosity'       , ...
                       'T'};
            
            model = dispatchParams(model, paramobj, fdnames);

            model.GasSupplyBc = ProtonicMembraneGasSupplyBc([]);
            
            model.constants = PhysicalConstants();
            
            model.gasInd.H2O = 1;
            model.gasInd.O2  = 2;
            model.nGas = 2;
            
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

        function state = updatePressure(model, state)

            nGas = model.nGas;
            
            p = state.pressures{1};

            for igas = 2 : nGas
                p = p + state.pressures{2};
            end

            state.pressure = p;
            
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
