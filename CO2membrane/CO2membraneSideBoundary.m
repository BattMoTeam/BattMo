classdef CO2membraneSideBoundary < BaseModel
    
    properties

        constants
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        
    end
    
    methods
        
        function model = CO2membraneSideBoundary(inputparams)
            
            model.constants = PhysicalConstants();
            
            model.gasInd.CO2 = 1;
            model.gasInd.O2  = 2;
            model.gasInd.N2  = 3;
            model.gasInd.Ar  = 4;
            model.nGas = 4;

        end

        
        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);

            nGas = model.nGas;
            
            varnames = {};
            %  fluxes [mol/s] (positive sign for outer fluxes)
            varnames{end + 1} = VarName({}, 'fluxes', nGas);
            %  mol fractions
            varnames{end + 1} = VarName({}, 'molFractions', nGas);
            %  total pressure [pascal]
            varnames{end + 1} = 'pressure';
            %  total flux
            varnames{end + 1} = 'flux';
            
            model = model.registerVarNames(varnames);

            fn = @CO2membraneSideBoundary.updateFluxes;
            inputvarnames = {VarName({}, 'molFractions', nGas), ...
                             'flux'};
            outputvarname = VarName({}, 'fluxes', nGas);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn = @CO2membraneSideBoundary.updateDependentMolFraction;
            inputvarnames = {VarName({}, 'molFractions', nGas, (1 : (nGas - 1))')};
            outputvarname = VarName({}, 'molFractions', nGas, nGas);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            
        end

        function state = updateDependentMolFraction(model, state)

            nGas = model.nGas;
            
            tmf = 0*state.molFractions{1};

            for igas = 1 : (nGas - 1)

                tmf = tmf + state.molFractions{igas};

            end

            state.molFractions{nGas} = 1 - tmf;

        end
        
        function state = updateFluxes(model, state)

            nGas = model.nGas;

            f = state.flux;
            
            for igas = 1 : nGas

                mf = state.molFractions{igas};
                fs{igas} = mf.*f;
            end

            state.fluxes = fs;
            
        end
        
    end
    
end

