classdef ProtonicMembraneGasSupplyBc < BaseModel
    
    properties

        molecularWeights
        T
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        constants
        
    end
    
    methods
        
        function model = ProtonicMembraneGasSupplyBc(inputparams)
            
            model = model@BaseModel();
            
            fdnames = {'molecularWeights', ...
                       'T'};

            model = dispatchParams(model, inputparams, fdnames);
            
            model.gasInd.H2O = 1;
            model.gasInd.O2  = 2;
            model.nGas = 2;

            model.constants = PhysicalConstants();
            
        end

        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);

            nGas = model.nGas;
            
            varnames = {};

            % Volume fractions
            varnames{end + 1} = VarName({}, 'massfractions', nGas);
            % Total pressure
            varnames{end + 1} = 'pressure';
            % Density
            varnames{end + 1} = 'density';
            % Densities
            varnames{end + 1} = VarName({}, 'densities', nGas);
            % Mass Flux terms at boundary (outward)
            varnames{end + 1} = VarName({}, 'massFluxes', nGas);
            % Boundary Equation relating boundary massFluxes and boundary pressure values
            varnames{end + 1} = VarName({}, 'bcFluxEquations', nGas);
            % Partial pressures
            varnames{end + 1} = VarName({}, 'pressures', nGas);
            
            model = model.registerVarNames(varnames);

            model = model.setAsExtraVarName(VarName({}, 'pressures', nGas));
            
            fn = @ProtonicMembraneGasSupply.updateMassFraction;
            inputvarnames = {VarName({}, 'massfractions', nGas, 1)};
            outputvarname = VarName({}, 'massfractions', nGas, 2);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
            fn = @(model, state) ProtonicMembraneGasSupply.updateDensity(model, state);
            fn = {fn, @(prop) PropFunction.literalFunctionCallSetupFn(prop)};
            inputvarnames = {VarName({}, 'massfractions', nGas), 'pressure'};
            model = model.registerPropFunction({'density', fn, inputvarnames});
            
            fn = @(model, state) ProtonicMembraneGasSupply.updateDensities(model, state);
            fn = {fn, @(prop) PropFunction.literalFunctionCallSetupFn(prop)};
            for igas = 1 : nGas
                inputvarnames = {'density', ...
                                 VarName({}, 'massfractions', nGas, igas)};
                outputvarname = VarName({}, 'densities', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            end

            fn = @(model, state) ProtonicMembraneGasSupply.updatePressures(model, state);
            fn = {fn, @(prop) PropFunction.literalFunctionCallSetupFn(prop)};
            inputvarnames = {VarName({}, 'massfractions', nGas), ...
                             'pressure'};
            outputvarname = VarName({}, 'pressures', nGas);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

        end
        
    end

    methods


        function state = updateMassFraction(model, state)

            state.massfractions{2} = 1 - state.massfractions{1};
            
        end

    end
    
end
