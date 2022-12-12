classdef PlatiniumCatalystLayer < CatalystLayer

    properties

    end
    
    methods
        
        function model = PlatiniumCatalystLayer(paramobj)

            model = model@CatalystLayer(paramobj);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            % Assemble the reaction rate constants
            fn = @() CatalystLayer.updateReactionRateConstants;
            inputnames = {'cOHelyte', 'H2OaElyte'};
            model = model.registerPropFunction({'elyteReactionRateConstant', fn, inputnames});
            inputnames = {'cOHinmr'};
            model = model.registerPropFunction({'inmrReactionRateConstant', fn, inputnames});
            
        end

        function state = updateReactionRateConstants(model, state)
            
            j0   = model.j0;
            cOH0 = model.elyteParams.OH.c0;
            
            cOH = state.cOHelyte;
            aw  = state.H2OaElyte;
            
            state.elyteReactionRateConstant = j0.*((cOH/cOH0.*aw).^0.5);
            
            aw = state.H2OaInmr;
            state.inmrReactionRateConstant = j0.*(aw.^0.5);
            
        end
        
    end
    
end
    
