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
            % ControlEquations
            varnames{end + 1} = 'controlEquations';
            
            model = model.registerVarNames(varnames);

            fn = @ProtonicMembraneGasSupply.updateVolumeFraction;
            inputvarnames = {VarName({}, 'volumefractions', nGas, 1)};
            outputvarname = VarName({}, 'volumefractions', nGas, 2);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
            fn = @ProtonicMembraneGasSupplyBc.updateControlEquations;
            inputvarnames = {'pressure'                          , ...
                             VarName({}, 'volumefractions', nGas), ...
                             VarName({}, 'massFluxes', nGas)};
            model = model.registerPropFunction({'controlEquations', fn, inputvarnames});
            
            fn = @(model, state) ProtonicMembraneGasSupply.updateDensity(model, state);
            fn = {fn, @(prop) PropFunction.literalFunctionCallSetupFn(prop)};
            inputvarnames = {VarName({}, 'volumefractions', nGas), 'pressure'};
            model = model.registerPropFunction({'density', fn, inputvarnames});
            
        end
        
    end

    methods

        function nbc = getNumberBcFaces(model)

            nbc = numel(model.controlHelpers.bcfaces);
            
        end


        function state = updateVolumeFraction(model, state)

            state.volumefractions{2} = 1 - state.volumefractions{1};
            
        end
        
        function state = updateControlEquations(model, state)

            controlHelpers = model.controlHelpers;

            maps   = controlHelpers.maps;
            bcvals = controlHelpers.vals;

            pressure   = state.pressure;
            vf         = state.volumefractions{1};
            massFluxes = state.massFluxes;
            
            eqs = {};
            for icomp = 1 : 2
                for itype = 1 : 2
                    switch itype
                      case 1
                        switch icomp
                          case 1
                            val = pressure;
                          case 2
                            val = vf;
                        end                            
                      case 2
                        val = massFluxes{icomp};
                      otherwise
                        error('type index not recognized');
                    end
                    eqs{end + 1} = maps{itype, icomp}*val - bcvals{itype, icomp};
                end
            end

            state.controlEquations = vertcat(eqs{:});
            
        end

    end
    
    methods(Static)

        function bceq = setupBcEquation(model, bcFlux, pIn, pBc, rhoIn, rhoBc, vfIn, vfBc, Tbc)
        % Tbc constains both transmissibility and (homogeneous) coefficient
            
            v  = Tbc.*(pIn - pBc);
            
            isOutward = (v > 0);

            v(isOutward)  = rhoIn(isOutward).*vfIn(isOutward).*v(isOutward);
            v(~isOutward) = rhoBc(~isOutward).*vfBc(~isOutward).*v(~isOutward);

            bceq = v - bcFlux;
            
        end
        
    end
    
end
