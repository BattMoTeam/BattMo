classdef ProtonicMembraneGasSupplyBc < BaseModel
    
    properties

        molecularWeights
        T
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        constants
        
        controlHelpers
        
    end
    
    methods
        
        function model = ProtonicMembraneGasSupplyBc(paramobj)
            
            model = model@BaseModel();
            
            model.gasInd.H2O = 1;
            model.gasInd.O2  = 2;
            model.nGas = 2;

            model.constants = PhysicalConstants();
            
        end

        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);

            nGas = model.nGas;
            
            varnames = {};

            % Gas partial pressures at boundary
            varnames{end + 1} = VarName({}, 'pressures', nGas);
            % Gas densities  at boundary
            varnames{end + 1} = VarName({}, 'densities', nGas);
            % Mass Flux terms at boundary (outward)
            varnames{end + 1} = VarName({}, 'massFluxes', nGas);
            % Boundary Equation relating boundary massFluxes and boundary pressure values
            varnames{end + 1} = VarName({}, 'boundaryEquations', nGas);
            % ControlEquations
            varnames{end + 1} = 'controlEquation';
            
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
                model = model.registerPropFunction({'controlEquation', fn, inputvarnames});
                
            end
            
        end
        
    end

    methods

        function nbc = getNumberBcFaces(model)

            nbc = numel(model.controlHelpers.bcfaces);
            
        end
            
        function state = updateControlEquations(model, state)

            controlHelpers = model.controlHelpers;

            maps = controlHelpers.maps;
            bcvals = controlHelpers.vals;

            pressures  = state.pressures;
            massFluxes = state.massFluxes;
            
            eqs = {};
            for icomp = 1 : 2
                for itype = 1 : 2
                    switch itype
                      case 1
                        val = pressures{icomp};
                      case 2
                        val = massFluxes{icomp};
                      otherwise
                        error('type index not recognized');
                    end
                    eqs{end + 1} = maps{itype, icomp}*val - bcvals{itype, icomp};
                end
            end

            state.controlEquation = vertcat(eqs{:});
            
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
