classdef MagnesiumActiveMaterial < SeaWaterActiveMaterial

    
    properties
        Asurf
    end
    
    methods

        function model = MagnesiumActiveMaterial(paramobj)
            
            model = model@SeaWaterActiveMaterial(paramobj);
            fdnames = {'Asurf', ...
                       'etaMax'};
            model = dispatchParams(model, paramobj, fdnames);
            
        end
        
        function model = registerVarAndPropfuncNames(model)
        
            model = registerVarAndPropfuncNames@SeaWaterActiveMaterial(model);
        
            varnames = {};
            
            varnames{end + 1} = 'specificSurfaceArea';
            
            model = model.registerVarNames(varnames);
            
            fn = @MagnesiumActiveMaterial.updateENernst;
            inputnames = {'T', 'cElectrolyte'};
            model = model.registerPropFunction({'ENernst', fn, inputnames});
            
            fn = @MagnesiumActiveMaterial.updateReactionRate;
            inputnames = {'eta', 'T', 'specificSurfaceArea'};
            model = model.registerPropFunction({'R', fn, inputnames}); 
            
        end
        
       
        function state = updateENernst(model, state)
            
            R = model.con.R;
            F = model.con.F;
            
            T = state.T;
            c = state.cElectrolyte;
            
            ENernst = -2.7 - R.*T./(2.*F).*log(1000./c);
            
            state.ENernst = ENernst;

        end
        
        function state = updateReactionRate(model, state)
            
            R = model.con.R;
            F = model.con.F;
            etaMax = model.etaMax;
            
            eta   = state.eta;
            T     = state.T;
            Asurf = state.specificSurfaceArea;
            
            R = Asurf./(2*F).*SmoothButlerVolmerEquation(1e-6, 0.5, 2, eta, T, etaMax); % mol m^-3 s^-1
            
            state.R = R;

        end            
        

    end
    
end

