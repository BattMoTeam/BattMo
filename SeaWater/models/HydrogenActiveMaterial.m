classdef HydrogenActiveMaterial < SeaWaterActiveMaterial

    properties
        Asurf % 
    end
    
    methods
        
        function model = HydrogenActiveMaterial(paramobj)
            
            model = model@SeaWaterActiveMaterial(paramobj);
            fdnames = {'Asurf'};
            model = dispatchParams(model, paramobj, fdnames);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@SeaWaterActiveMaterial(model);
            
            fn = @() HydrogenActiveMaterial.updateENernst;
            inputnames = {'T', 'cElectrolyte'};
            model = model.registerPropFunction({'ENernst', fn, inputnames});
            
            fn = @() HydrogenActiveMaterial.updateReactionRate;
            inputnames = {'eta', 'T'};
            model = model.registerPropFunction({'R', fn, inputnames});
            
        end
        
        function state = updateENernst(model, state)
            
            R = model.con.R;
            F = model.con.F;
            
            T = state.T;
            c = state.cElectrolyte;
            
            ENernst = 0 - R.*T./(2*F).*log(1000^2./c.^2);
            
            state.ENernst = ENernst;
        end
                
        function state = updateReactionRate(model, state)

            F     = model.con.F;
            Asurf = model.Asurf;
            
            eta = state.eta;
            T = state.T;
            
            R = Asurf./(2.*F).*SmoothButlerVolmerEquation(5e-4, 0.5, 2, eta, T, inf);
            
            state.R = R;
        end
                
    end
    
    
end
