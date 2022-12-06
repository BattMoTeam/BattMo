classdef Catalyser < BaseModel

    properties

        constants

        compInd   % mapping structure for component indices
        phaseInd  % mapping structure for phase indices
        liquidInd % mapping structure for component indices
        gasInd    % mapping structure for component indices

        E0 % Standard equilibrium potential,   [V]
        j0 % Exchange current density
        inmr %
        % inmr.cT   % Total concentration of charged groups
        % inmr.kxch % Rate constant for exchange between ionomer and electrolyte. [s^-1]
        % inmr.OH.z
        PF1 % should be initialized using zero

    end

    methods

        function model = Catalyser(paramobj)

        end

        function model = registerVarAndPropfuncNames(model)

        %% Declaration of the Dynamical Variables and Function of the model
        % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            names = {};
            % Temperature
            varnames = {'T'};
            % Catalyser electrical potential
            varnames{end + 1} = {'phi'};
            % Electric Potential in electrolyte and ionomer
            varnames{end + 1} = 'phiElyte';
            varnames{end + 1} = 'phiInmr';
            % concentration of OH in electrotyte and ionomer
            varnames{end + 1} = 'cOHelyte';
            varnames{end + 1} = 'cOHinmr';
            % Equilibrium Potential for electrolyte and ionomer
            varnames{end + 1} = 'Eelyte';
            varnames{end + 1} = 'Einmr';
            % Partial pressure of the active gas (H2 or O2 for exampler) in the electrolyte
            varnames{end + 1} = 'pressureActiveGas';
            % Potential difference with respect to the electrolyte and ionomer
            varnames{end + 1} = {'etaElyte'};
            varnames{end + 1} = {'etaInmr'};
            % Water activity in electrolyte and ionomer
            varnames{end + 1} = 'H2OaElyte';
            varnames{end + 1} = 'H2OaInmr';
            % Reaction rate constant (j0) for electrolyte and ionomer
            varnames{end + 1} = 'elyteReactionRateConstant';
            varnames{end + 1} = 'inmrReactionRateConstant';
            % Reaction rate for electrolyte and ionomer
            varnames{end + 1} = 'elyteReactionRate';
            varnames{end + 1} = 'inmrReactionRate';

            model = model.registerVarNames(varnames);

            % Assemble equilibrium Potential for electrolyte
            fn = @() Catalyser.updateEelyte;
            inputnames = {'T', 'cOHelyte', 'pressureActiveGas', 'H2OaElyte'};
            model = model.registerPropFunction({'Eelyte', fn, inputnames});

            % Assemble equilibrium Potential for inmr
            fn = @() Catalyser.updateEinmr;
            inputnames = {'T', 'cOHinmr', 'pressureActiveGas', 'H2OaInmr'};
            model = model.registerPropFunction({'Einmr', fn, inputnames});

            % Assemble reactive potential
            fn = @() Catalyser.updateEtas;
            inputnames = {'phiElyte', 'phi' 'Eelyte'};
            model = model.registerPropFunction({'etaElyte', fn, inputnames});
            inputnames = {'phiInmr', 'phi' 'Einmr'};
            model = model.registerPropFunction({'etaInmr', fn, inputnames});            

            % Assemble the reaction rate constants
            fn = @() Catalyser.updateReactionRateConstants;
            inputnames = {'cOHelyte'};
            model = model.registerPropFunction({'elyteReactionRateConstant', fn, inputnames});
            inputnames = {'cOHinmr'};
            model = model.registerPropFunction({'inmrReactionRateConstant', fn, inputnames});
            
            % Assemble the reaction rates
            fn = @() Catalyser.updateReactionRates;
            inputnames = {'elyteReactionRateConstant', 'etaElyte'};
            model = model.registerPropFunction({'elyteReactionRate', fn, inputnames});
            inputnames = {'inmrReactionRateConstant', 'etaInmr'};
            model = model.registerPropFunction({'inmrReactionRate', fn, inputnames});            

        end

        function state = updateEr(model, state)

            T = state.T;
            cHElyte = state.CHElyte;

            E0  = model.E0;
            con = model.constants;
            F   = con.F;
            c0  = con.c0;
            R   = con.R;

            state.Er = E0 - R.*T./(2.*F) .* log(c0.^2 .* cHElyte.^-2);
        end

        function state = updateSorption(model, state)

            H2OaElyte = state.H2OaElyte;
            H2OaInmr = state.H2OaInmr;

            kML = model.inmr.kML;

            state.Rsorption = ionomerSorption(kML, H2OaElyte, H2OaInmr);
        end

        function state = updateEtas(model, state)

            phiElyte = state.phiElyte;
            phiInmr  = state.phiInmr;
            E        = state.E;
            Er       = state.Er;

            state.eta      = phiInmr - phiElyte;
            state.etaElyte = E - phiElyte - Er;
            state.etaInmr  = E - phiInmr - Er;
        end

        function state = updateIonomerExchange(model, state)

            eta      = state.eta;
            T        = state.T;
            cOHElyte = state.cOHElyte;

            ir = model.inmr;

            state.Rxch = ionomerExchange(ir.kxch, ir.cT, ir.OH.z, eta, T, cOHElyte);
        end

        function state = updateFluxes(model, state)

         etaElyte = state.etaElyte;
         etaInmr  = state.etaInmr;
         H2OaInmr = state.H2OaInmr;
         T        = state.T;
         gasp1    = state.gasPressuresElyte{1};
         gasp2    = state.gasPressuresElyte{2};

         ir    = model.inmr;
         Asp   = model.Asp;
         PF1   = model.PF1;
         j0    = model.j0;
         alpha = model.alpha;
         n     = model.n;

         X_inmr = 0.99;

         PF = (ir.cT/1000).*(gasp1./1e5).*(gasp2./1e5);
         jElyte = (1 - X_inmr).*Asp'.*butlerVolmer(j0, alpha, n, etaElyte, T);

         PF2 = H2OaInmr;
         jInmr = X_inmr.*Asp'*butlerVolmer_split(j0, PF1, PF2, alpha, n, etaInmr, T);

         state.jElyte = jElyte;
         state.jInmr  = jInmr;
         state.j      = jElyte + jInmr;
        end


    end

end
