classdef PlatiniumCatalystLayer < CatalystLayer

    properties

    end
    
    methods
        
        function model = PlatiniumCatalystLayer(paramobj)

            model = model@CatalystLayer(paramobj);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@CatalystLayer(model);

            % Assemble the reaction rate constants
            fn = @() PlatiniumCatalystLayer.updateReactionRateConstants;
            inputnames = {'cOHelyte', 'H2OaElyte'};
            model = model.registerPropFunction({'elyteReactionRateConstant', fn, inputnames});
            inputnames = {'cOHinmr'};
            model = model.registerPropFunction({'inmrReactionRateConstant', fn, inputnames});
            

            % Assemble equilibrium Potential for electrolyte
            fn = @() PlatiniumCatalystLayer.updateEelyte;
            inputnames = {'T', 'cOHelyte', 'pressureActiveGas', 'H2OaElyte'};
            model = model.registerPropFunction({'Eelyte', fn, inputnames});

            % Assemble equilibrium Potential for inmr
            fn = @() PlatiniumCatalystLayer.updateEinmr;
            inputnames = {'T', 'cOHinmr', 'pressureActiveGas', 'H2OaInmr'};
            model = model.registerPropFunction({'Einmr', fn, inputnames});

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
        
        function state =  updateSources(model, state)

        % Reaction in Electrolyte : 2*H2O + 2*e- <-> + H2 + 2(OH-)_elyte 
        % Reaction in Membrane :    2*H2O + 2*e- <-> + H2 + 2(OH-)_inmr

            F = model.constants.F;
            vols = model.G.cells.volumes;
            
            elyteR = state.elyteReactionRate;
            inmrR  = state.inmrReactionRate;

            R = elyteR + inmrR;
            
            state.activeGasSource = R/F;
            state.elyteH2Osource  = -2*R/F;
            state.elyteOHSource   = 2*elyteR/F;
            state.inmrOHSource    = 2*inmrR/F;
            state.eSource         = -2*R.*vols;
           
        end

         function state = updateEelyte(model, state)

            T    = state.T;
            cOH  = state.cOHElyte;
            H2Oa = state.H2OaElyte;
            
            E0  = model.E0;
            c0  = model.sp.OH.c0;
            con = model.constants;

            F  = con.F;
            R  = con.R;

            state.Eelyte = E0 + R*T/(2*F)*log((H2Oa.^2).*(c0.^2).*(cOH.^-2));

         end
 
         function state = updateEinmr(model, state)

            T    = state.T;
            cOH  = state.cOHinmr;
            H2Oa = state.H2OaInmr;
            
            E0  = model.E0;
            c0  = model.sp.OH.c0;
            con = model.constants;

            F  = con.F;
            R  = con.R;

            state.Einmr = E0 + R*T/(2*F)*log((H2Oa.^2).*(c0.^2).*(cOH.^-2));
            
         end
        
    end
    
end
    
