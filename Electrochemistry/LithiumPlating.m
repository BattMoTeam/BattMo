classdef LithiumPlating < BaseModel

    properties
        alphaA1 = 0.5
        alphaC1 = 0.5 % Transfer coefficients of lithium intercalation reaction
        T = 298.15           % Temperature [K]
        U = 0;               % Open Circuit Voltage (OCV)
        k2 = 1e-10;         % Reaction rate constant [A m^-2 (mol m^-3)^-1.5]
        cmax = 1000;         % Maximum lithium concentration [mol/m³]
        Uref % depends on position

        beta = 1000 %need to check again
        reversibleFraction = 0.775
        deadFraction = 0.175
        SEIFraction = 0.05
    
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
                        'platingFlux', ...
                        'strippingFlux', ...
                        'massAccum'     , ...
                        'massCons'
                        'U'
                        'exchangeCurrentDensity'
                        'overpotential'
                        'nLiRev'
                        };

            model = model.registerVarNames(varnames);
            
            fn = @LithiumPlating.updateMassCons; 
            model = model.registerPropFunction({'massCons', fn, {'massAccum', 'flux'}});
            
            fn = @LithiumPlating.updateMassAccum; 
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({'massAccum', fn, {'concentration'}});

            fn = @LithiumPlating.updateExchangeCurrentDensity; 
            model = model.registerPropFunction({'exchangeCurrentDensity', fn, {'cElectrolyte', 'concentration'}});
            
            fn = @LithiumPlating.updateOverpotential; 
            model = model.registerPropFunction({'overpotential', fn, {'phiElectrode', 'phiElectrolyte'}});
     
            fn = @LithiumPlating.updatePlatingFlux; 
            model = model.registerPropFunction({'platingFlux', fn, {'overpotential', 'exchangeCurrentDensity'}});
            
            fn = @LithiumPlating.updateStrippingFlux; 
            model = model.registerPropFunction({'strippingFlux', fn, {'overpotential', 'exchangeCurrentDensity'}});

            fn = @LithiumPlating.updateNLiRev; 
            model = model.registerPropFunction({'nLiRev', fn, {'platingFlux'}});
            
        end
        
        function state = updateExchangeCurrentDensity(model, state) %eq 20
            
            k2 = model.k2
            ce = state.cElectrolyte;
            alphaA2 = model.alphaA2
            
            i0 = k2*ce^alphaA2

            state.exchangeCurrentDensity = i0

        end

        function state = updateOverpotential(model, state) %must be completed
            eta = state.phiElectrode - state.phiElectrolyte
            state.overpotential = eta
        end

        function state = updatePlatingFlux(model, state)
            eta = state.overpotential
            if eta > 0 %checking if lithium is plating
                state.platingFlux = 0;
            else
                F = model.constants.F;
                R = model.constants.R;
                T = model.T;
                alphaA2 = model.alphaA2
                alphaC2 = model.alphaC2
                i0 = state.exchangeCurrentDensity

                %butler-Volmer equation (current density, A/m²)
                j = i0 .* (exp((alphaA2 * F * eta) / (R * T)) - exp((-alphaC2 * F * eta) / (R * T)));

                %convert to molar flux (mol/m²/s) by dividing by Faraday constant
                platingFlux = j ./ F;

                state.platingFlux = platingFlux;
            end
        end

        function state = updateStrippingFlux(model, state)
            eta = state.overpotential
            if eta < 0 | state.nLiRev <= 0 %checking if lithium is stripping
                state.platingFlux = 0;
            else
                F = model.constants.F;
                R = model.constants.R;
                T = model.T;
                alphaA2 = model.alphaA2
                alphaC2 = model.alphaC2
                i0 = state.exchangeCurrentDensity

                %correction factor considering the limitation by the amount of reversible lithium
                
                correction = model.beta.*state.nLiRev./(1 + model.beta.*state.nLiRev)

                %butler-Volmer equation (current density, A/m²)
                j = correction .* i0 .* (exp((alphaA2 * F * eta) / (R * T)) - exp((-alphaC2 * F * eta) / (R * T)));

                %convert to molar flux (mol/m²/s) by dividing by Faraday constant
                platingFlux = j ./ F;

                state.platingFlux = platingFlux;
            end
        end

        function state = updateMassAccum(model, state, state0, dt)
            state.massAccum = (state.concentration - state0.concentration) ./ dt;         
        end

        function state = updateMassCons(model, state) %need to see that again
            flux = state.flux;
            accum = state.massAccum;
            cons = assembleConservationEquation(model, flux, 0, 0, accum);
            state.massCons = cons;
        end
    
    end

end
