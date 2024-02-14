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
            varnames{end + 1} = VarName({}, 'volumefractions', nGas);
            % Total pressure
            varnames{end + 1} = 'pressure';
            % Density
            varnames{end + 1} = 'density';
            % Mass Flux terms at boundary (outward)
            varnames{end + 1} = VarName({}, 'massFluxes', nGas);
            % Boundary Equation relating boundary massFluxes and boundary pressure values
            varnames{end + 1} = VarName({}, 'boundaryEquations', nGas);
            
            model = model.registerVarNames(varnames);

            fn = @ProtonicMembraneGasSupply.updateVolumeFraction;
            inputvarnames = {VarName({}, 'volumefractions', nGas, 1)};
            outputvarname = VarName({}, 'volumefractions', nGas, 2);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
            
            fn = @(model, state) ProtonicMembraneGasSupply.updateDensity(model, state);
            fn = {fn, @(prop) PropFunction.literalFunctionCallSetupFn(prop)};
            inputvarnames = {VarName({}, 'volumefractions', nGas), 'pressure'};
            model = model.registerPropFunction({'density', fn, inputvarnames});
            
        end
        
    end

    methods


        function state = updateVolumeFraction(model, state)

            state.volumefractions{2} = 1 - state.volumefractions{1};
            
        end

    end
    
end
