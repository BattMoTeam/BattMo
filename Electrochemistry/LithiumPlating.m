classdef LithiumPlating < BaseModel

    properties
        alpha = 0.5          % Symmetry factor for Butler-Volmer equation
        T = 298.15           % Temperature [K]
        U = 0;               % Open Circuit Voltage (OCV)
        kLi = 1e-10;         % Reaction rate constant [A m^-2 (mol m^-3)^-1.5]
        cmax = 1000;         % Maximum lithium concentration [mol/m³]
    end

    methods

        function model = LithiumPlating(inputparams)

            model = model@BaseModel();
            
            fdnames = {'alpha', ...
                       'kLi', ...
                       'U', ...
                       'cmax'};
            model = dispatchParams(model, inputparams, fdnames);
            
        end

        function model = registerVarAndPropfuncNames(model)

            varnames = {'phiElectrode'  , ...
                        'phiElectrolyte', ...
                        'cElectrolyte'  , ...
                        'concentration' , ...
                        'flux'          , ...
                        'massAccum'     , ...
                        'massCons'};

            model = model.registerVarNames(varnames);
            
            fn = @LithiumPlating.updateMassCons; 
            model = model.registerPropFunction({'massCons', fn, {'massAccum', 'flux'}});
            
            fn = @LithiumPlating.updateMassAccum; 
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({'massAccum', fn, {'concentration'}});

            fn = @LithiumPlating.updateFlux; 
            model = model.registerPropFunction({'flux', fn, {'phiElectrode', 'phiElectrolyte', 'cElectrolyte', 'concentration'}});
        end

        function state = updateFlux(model, state)
            
            F = model.constants.F; 
            R = model.constants.R;
            T = model.T;

            phi_s = state.phiElectrode;
            phi_e = state.phiElectrolyte;

            ce = state.cElectrolyte;
            cs = state.concentration;

            U = model.U;

            eta = phi_s - phi_e - U; % Overpotential

            k0 = model.kLi;
            alpha = model.alpha;
            cmax = model.cmax;

            %exchange current density j0
            j0 = k0 .* sqrt(ce) .* sqrt(cs) .* sqrt(cmax - cs);

            %butler-Volmer equation (current density, A/m²)
            j = j0 .* (exp((alpha * F * eta) / (R * T)) - exp((-(1 - alpha) * F * eta) / (R * T)));

            %convert to molar flux (mol/m²/s) by dividing by Faraday constant
            flux = j ./ F;

            state.flux = flux;
        end

        function state = updateMassAccum(model, state, state0, dt)
            state.massAccum = (state.concentration - state0.concentration) ./ dt;         
        end

        function state = updateMassCons(model, state)
            flux = state.flux;
            accum = state.massAccum;
            cons = assembleConservationEquation(model, flux, 0, 0, accum);
            state.massCons = cons;
        end
    
    end

end
