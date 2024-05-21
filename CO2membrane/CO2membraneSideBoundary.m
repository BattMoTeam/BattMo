classdef CO2membraneSideBoundary < BaseModel
    
    properties

        constants
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        
    end
    
    methods
        
        function model = CO2membraneSideBoundary(inputparams)
            
            model.constants = PhysicalConstants();

            model = CO2membrane.setupGasStructures(model);

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


            for igas = 1 : nGas
                fn = @CO2membraneSideBoundary.updateFluxes;
                inputvarnames = {VarName({}, 'molFractions', nGas, igas), ...
                                 'flux'};
                outputvarname = VarName({}, 'fluxes', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            end
            
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

