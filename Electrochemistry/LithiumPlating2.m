classdef LithiumPlatingLatz < BaseModel

    properties
        T = 298.15
        F = 96485
        R = 8.314

        alphaPl = 0.3
        alphaStr = 0.7
        alphaChInt = 0.5

        kPl = 1e-10
        kChInt = 1e-12

        nPl0 = 1e-6
        muLiRef = 0

        nPlLimit = 2.3e-5

        SEIFraction = 0.05
        MSEI = 0.162
        rhoSEI = 1690
        deltaSEI0 = 1e-9
        sigmaSEI = 5e-6
    end

    methods

        function model = LithiumPlatingLatz(inputparams)
            model = model@BaseModel();
            fdnames = {};
            model = dispatchParams(model, inputparams, fdnames);
        end

        function model = registerVarAndPropfuncNames(model)

            varnames = {'phiSolid'        , ...
                        'phiElectrolyte'  , ...
                        'cElectrolyte'    , ...
                        'nPl'             , ...
                        'nPlAccum'        , ...
                        'nPlCons'         , ...
                        'platingFlux'     , ...
                        'chemicalFlux'    , ...
                        'etaPlating'      , ...
                        'etaChemical'     , ...
                        'activityPlated'  , ...
                        'surfaceCoverage' , ...
                        'nSEI'            , ...
                        'nSEIAccum'       , ...
                        'nSEICons'        , ...
                        'SEIThickness'    };

            model = model.registerVarNames(varnames);

            fn = @LithiumPlatingLatz.updateNPlAccum;
            fn = {fn, @(pf) PropFunction.accumFuncCallSetupFn(pf)};
            model = model.registerPropFunction({'nPlAccum', fn, {'nPl'}});

            fn = @LithiumPlatingLatz.updateActivityPlated;
            model = model.registerPropFunction({'activityPlated', fn, {'nPl'}});

            fn = @LithiumPlatingLatz.updateEtaPlating;
            model = model.registerPropFunction({'etaPlating', fn, {'phiSolid', 'phiElectrolyte', 'activityPlated', 'SEIThickness'}});

            fn = @LithiumPlatingLatz.updateEtaChemical;
            model = model.registerPropFunction({'etaChemical', fn, {'activityPlated'}});

            fn = @LithiumPlatingLatz.updatePlatingFlux;
            model = model.registerPropFunction({'platingFlux', fn, {'cElectrolyte', 'etaPlating'}});

            fn = @LithiumPlatingLatz.updateChemicalFlux;
            model = model.registerPropFunction({'chemicalFlux', fn, {'etaChemical'}});

            fn = @LithiumPlatingLatz.updateNPlCons;
            model = model.registerPropFunction({'nPlCons', fn, {'nPlAccum', 'platingFlux', 'chemicalFlux'}});

            fn = @LithiumPlatingLatz.updateSurfaceCoverage;
            model = model.registerPropFunction({'surfaceCoverage', fn, {'nPl'}});

            fn = @LithiumPlatingLatz.updateNSEIAccum;
            fn = {fn, @(pf) PropFunction.accumFuncCallSetupFn(pf)};
            model = model.registerPropFunction({'nSEIAccum', fn, {'nSEI'}});

            fn = @LithiumPlatingLatz.updateNSEICons;
            model = model.registerPropFunction({'nSEICons', fn, {'platingFlux', 'nSEIAccum'}});

            fn = @LithiumPlatingLatz.updateSEIThickness;
            model = model.registerPropFunction({'SEIThickness', fn, {'nSEI'}});
        end

        function state = updateNPlAccum(model, state, state0, dt)
            state.nPlAccum = (state.nPl - state0.nPl) / dt;
        end

        function state = updateActivityPlated(model, state)
            nPl = state.nPl;
            n0 = model.nPl0;
            state.activityPlated = nPl ./ (nPl + n0);
        end

        function state = updateEtaPlating(model, state)
            phiS = state.phiSolid;
            phiE = state.phiElectrolyte;
            aPl = state.activityPlated;

            RSEI = state.SEIThickness / model.sigmaSEI;
            j = model.kPl*0;  % !!! Approximation: use i0 as constant for RSEI drop

            eta = phiS - phiE - j * model.F * RSEI + (model.R * model.T / model.F) * log(aPl);
            state.etaPlating = eta;
        end

        function state = updateEtaChemical(model, state)
            aPl = state.activityPlated;
            eta = -(model.R * model.T / model.F) * log(aPl);
            state.etaChemical = eta;
        end

        function state = updatePlatingFlux(model, state)
            eta = state.etaPlating;
            ce = state.cElectrolyte;
            T = model.T; R = model.R; F = model.F;

            i0 = model.kPl * ce.^model.alphaPl;
            j = i0 .* (exp((model.alphaPl * F * eta) / (R * T)) - ...
                       exp((-model.alphaStr * F * eta) / (R * T)));

            state.platingFlux = j ./ F;
        end

        function state = updateChemicalFlux(model, state)
            eta = state.etaChemical;
            jCh = model.kChInt * (exp(0.5 * model.F * eta / (model.R * model.T)) - ...
                                  exp(-0.5 * model.F * eta / (model.R * model.T)));
            state.chemicalFlux = jCh ./ model.F;
        end

        function state = updateNPlCons(model, state)
            flux = state.platingFlux - state.chemicalFlux;
            accum = state.nPlAccum;
            state.nPlCons = assembleConservationEquation(model, flux, 0, 0, accum);
        end

        function state = updateSurfaceCoverage(model, state)
            nPl = state.nPl;
            state.surfaceCoverage = min(nPl ./ model.nPlLimit, 1.0);
        end

        function state = updateNSEIAccum(model, state, state0, dt)
            state.nSEIAccum = (state.nSEI - state0.nSEI) / dt;
        end

        function state = updateNSEICons(model, state)
            flux = state.platingFlux * model.SEIFraction;
            accum = state.nSEIAccum;
            state.nSEICons = assembleConservationEquation(model, flux, 0, 0, accum);
        end

        function state = updateSEIThickness(model, state)
            delta = model.deltaSEI0 + (model.MSEI * state.nSEI) / (model.rhoSEI * model.F);
            state.SEIThickness = delta;
        end

    end
end
