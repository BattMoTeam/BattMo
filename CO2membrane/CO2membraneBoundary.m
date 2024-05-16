classdef CO2membraneBoundary < BaseModel
    
    properties

        Boundary
        
        constants
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        
    end
    
    methods
        
        function model = CO2membraneBoundary(inputparams)
            
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
            % feed fluxes [mol/s]
            varnames{end + 1} = VarName({}, 'feedFluxes', nGas);
            % permeate fluxes [mol/s]
            varnames{end + 1} = VarName({}, 'permeateFluxes', nGas);
            % feed pressures [pascal]
            varnames{end + 1} = VarName({}, 'feedPressures', nGas);
            % permeate pressures [pascal]
            varnames{end + 1} = VarName({}, 'permeatePressures', nGas);
            % feed total pressure [pascal]
            varnames{end + 1} = 'feedPressure';
            % permeate total pressure [pascal]
            varnames{end + 1} = 'permeatePressure';
            
            model = model.registerVarNames(varnames);

            fn = @CO2membrane.updateFeedPressure;
            inputvarnames = {VarName({}, 'feedPressures', nGas)};
            outputvarname = 'feedPressure';
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
            fn = @CO2membrane.updatePermeatePressure;
            inputvarnames = {VarName({}, 'permeatePressures', nGas)};
            outputvarname = 'permeatePressure';
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
        end

        
        function state = updateFeedPressure(model, state)

            nGas = model.nGas;
            
            p = 0*state.feedPressures{1};

            for igas = 1 : nGas
                p = p + state.feedPressures{igas};
            end

            state.feedPressure = p;
            
        end
        
        function state = updatePermeatePressure(model, state)

            nGas = model.nGas;
            
            p = 0*state.permeatePressures{1};

            for igas = 1 : nGas
                p = p + state.permeatePressures{igas};
            end

            state.permeatePressure = p;
            
        end
    end
    
end

