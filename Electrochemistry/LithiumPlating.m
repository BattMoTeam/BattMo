classdef LithiumPlatingLatz < BaseModel

    properties
        T = 298.15            % Temperature [K]
        F = 96485             % Faraday constant [C/mol]
        R = 8.314             % Universal gas constant [J/mol/K]

        alphaPl = 0.3         % Anodic transfer coefficient for plating/stripping
        alphaStr = 0.7        % Cathodic transfer coefficient for plating/stripping
        alphaChInt = 0.5      % Symmetry factor for chemical intercalation

        kPl = 1e-10           % Kinetic rate constant for plating [mol/(m²·s·(mol/m³)^α)]
        kChInt = 1e-12        % Kinetic rate constant for chemical intercalation [mol/(m²·s)]

        nPl0 = 1e-6           % Reference plated lithium amount for activity expression [mol/m²]
        muLiRef = 0           % Reference chemical potential of lithium metal [J/mol]

        nPlLimit = 2.3e-5     % Maximum plated Li before full coverage (1 monolayer) [mol/m²]
        ATotal = 1            % Total available surface area per unit volume [m²]
    end

    methods

        function model = LithiumPlatingLatz(inputparams)
            model = model@BaseModel();
            fdnames = {};
            model = dispatchParams(model, inputparams, fdnames);
        end

        function model = registerVarAndPropfuncNames(model)

            varnames = {'phiSolid'         , ...
                        'phiElectrolyte'   , ...
                        'cElectrolyte'     , ...
                        'nPl'              , ...
                        'nPlAccum'         , ...
                        'nPlCons'          , ...
                        'platingFlux'      , ...
                        'chemicalFlux'     , ...
                        'etaPlating'       , ...
                        'etaChemical'      , ...
                        'activityPlated'   , ...
                        'surfaceCoverage'  };

            model = model.registerVarNames(varnames);

            fn = @LithiumPlatingLatz.updateNPlAccum;
            fn = {fn, @(pf) PropFunction.accumFuncCallSetupFn(pf)};
            model = model.registerPropFunction({'nPlAccum', fn, {'nPl'}});

            fn = @LithiumPlatingLatz.updateActivityPlated;
            model = model.registerPropFunction({'activityPlated', fn, {'nPl'}});

            fn = @LithiumPlatingLatz.updateEtaPlating;
            model = model.registerPropFunction({'etaPlating', fn, {'phiSolid', 'phiElectrolyte', 'activityPlated'}});

            fn = @LithiumPlatingLatz.updateEtaChemical;
            model = model.registerPropFunction({'etaChemical', fn, {'activityPlated'}});

            fn = @LithiumPlatingLatz.updatePlatingFlux;
            model = model.registerPropFunction({'platingFlux', fn, {'cElectrolyte', 'etaPlating'}});

            fn = @LithiumPlatingLatz.updateChemicalFlux;
            model = model.registerPropFunction({'chemicalFlux', fn, {'etaChemical'}});

            fn = @LithiumPlatingLatz.updateNPlCons;
            model = model.registerPropFunction({'nPlCons', fn, {'nPlAccum', 'platingFlux', 'chemicalFlux', 'SurfaceCoverage'}});

            fn = @LithiumPlatingLatz.updateSurfaceCoverage;
            model = model.registerPropFunction({'surfaceCoverage', fn, {'nPl'}});
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

            eta = phiS - phiE + (model.R * model.T / model.F) * log(aPl);
            state.etaPlating = eta;
        end

        function state = updateEtaChemical(model, state)
            % Equation (23) of Hein et al.
            aPl = state.activityPlated;
            eta = -(model.R * model.T / model.F) * log(aPl); % Assuming U0 = 0
            state.etaChemical = eta;
        end

        function state = updatePlatingFlux(model, state)
            eta = state.etaPlating;
            ce = state.cElectrolyte;
            T = model.T; R = model.R; F = model.F;

            i0 = model.kPl * ce.^model.alphaPl;
            j = i0 .* (exp((model.alphaPl * F * eta) / (R * T)) - ...
                       exp((-model.alphaStr * F * eta) / (R * T)));

            state.platingFlux = j ./ F;  % [mol/m²/s]
        end

        function state = updateChemicalFlux(model, state)
            eta = state.etaChemical;
            jCh = model.kChInt * (exp(0.5 * model.F * eta / (model.R * model.T)) - ...
                                  exp(-0.5 * model.F * eta / (model.R * model.T)));
            state.chemicalFlux = jCh ./ model.F; % [mol/m²/s]
        end

        function state = updateNPlCons(model, state)
            flux = state.platingFlux- state.chemicalFlux;
            accum = state.nPlAccum;
            cons = assembleConservationEquation(model, flux, 0, 0, accum);
            state.nPlCons = cons;
        end

        function state = updateSurfaceCoverage(model, state)
            nPl = state.nPl;
            state.surfaceCoverage = min(nPl ./ model.nPlLimit, 1.0);
        end

    end
end
