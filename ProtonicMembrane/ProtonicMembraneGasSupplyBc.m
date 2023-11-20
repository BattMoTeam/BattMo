classdef ProtonicMembraneGasSupplyBc < BaseModel
    
    properties

        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.
        
    end
    
    methods
        
        function model = ProtonicMembraneGasSupplyBc(paramobj)
            
            model = model@BaseModel();
            
            model.gasInd.H2O = 1;
            model.gasInd.O2  = 2;
            model.nGas = 2;
            
        end

        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);

            nGas = model.nGas;
            
            varnames = {};

            % Gas partial pressures at boundary
            varnames{end + 1} = VarName({}, 'pressures', nGas);
            % Gas densities  at boundary
            varnames{end + 1} = VarName({}, 'densities', nGas);
            % Mass Flux terms at boundary
            varnames{end + 1} = VarName({}, 'massFluxes', nGas);
            % Boundary Equation relating boundary massFluxes and boundary pressure values
            varnames{end + 1} = VarName({}, 'boundaryEquations', nGas);
            % ControlEquations
            varnames{end + 1} = VarName({}, 'controlEquations', nGas);
            
            model = model.registerVarNames(varnames);

            for igas = 1 : nGas
                
                fn = @(model, state) ProtonicMembraneGasSupply.updateDensities(model, state);
                fn = {fn, @(prop) PropFunction.literalFunctionCallSetupFn(prop)};
                inputvarnames = {VarName({}, 'pressures', nGas, igas)};
                outputvarname = VarName({}, 'densities', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @ProtonicMembraneGasSupplyBc.updateControlEquations;
                inputvarnames = {VarName({}, 'pressures', nGas, igas), ...
                                 VarName({}, 'massFluxes', nGas, igas)};
                outputvarname = VarName({}, 'controlEquations', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
            end
            
        end
        
    end
    
    methods(Static)

        function bceq = setupBcEquation(model, bcFlux, pIn, pBc, rhoIn, rhoBc, Tbc)
        % Tbc constains both transmissibility and (homogeneous) coefficient
            
            dp = pIn - pBc;
            v  = Tbc.*dp;
            
            isOutward = (v > 0);

            v(isOutward)  = rhoIn(isOutward).*v(isOutward);
            v(~isOutward) = rhoBc(~isOutward).*v(~isOutward);

            bceq = v - bcFlux;
            
        end
        
    end
    
end
