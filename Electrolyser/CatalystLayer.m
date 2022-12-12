classdef CatalystLayer < BaseModel

    properties

        constants

        compInd   % mapping structure for component indices
        phaseInd  % mapping structure for phase indices
        liquidInd % mapping structure for component indices
        gasInd    % mapping structure for component indices

        j0 % Exchange current density
        inmrParams 
        % inmrParams.OH.z  : Charge
        % inmrParams.OH.c0 : OH reference concentration
        % inmrParams.E0    : standard equilibrium potential
        elyteParams
        % elyteParams.OH.c0 : reference OH concentration
        % elyteParams.E0    : standard equilibrium potential
        
        Xinmr % Fraction of specific area that is coversed with ionomer [-]
        volumetricSurfaceArea % Volumetric surface area [m^ -1]
        
    end

    methods

        function model = CatalystLayer(paramobj)

            model = model@BaseModel();
            
        end

        function model = registerVarAndPropfuncNames(model)

        %% Declaration of the Dynamical Variables and Function of the model
        % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            names = {};
            % Temperature
            varnames = {'T'};
            % CatalystLayer electrical potential
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
            fn = @() CatalystLayer.updateEelyte;
            inputnames = {'T', 'cOHelyte', 'pressureActiveGas', 'H2OaElyte'};
            model = model.registerPropFunction({'Eelyte', fn, inputnames});

            % Assemble equilibrium Potential for inmr
            fn = @() CatalystLayer.updateEinmr;
            inputnames = {'T', 'cOHinmr', 'pressureActiveGas', 'H2OaInmr'};
            model = model.registerPropFunction({'Einmr', fn, inputnames});

            % Assemble reactive potential
            fn = @() CatalystLayer.updateEtas;
            inputnames = {'phiElyte', 'phi' 'Eelyte'};
            model = model.registerPropFunction({'etaElyte', fn, inputnames});
            inputnames = {'phiInmr', 'phi' 'Einmr'};
            model = model.registerPropFunction({'etaInmr', fn, inputnames});            


            % Assemble the reaction rates
            fn = @() CatalystLayer.updateReactionRates;
            inputnames = {'elyteReactionRateConstant', 'etaElyte'};
            model = model.registerPropFunction({'elyteReactionRate', fn, inputnames});
            inputnames = {'inmrReactionRateConstant', 'etaInmr'};
            model = model.registerPropFunction({'inmrReactionRate', fn, inputnames});            

        end

        function state = updateEelyte(model, state)

            T    = state.T;
            cOH  = state.cOHelyte;
            H2Oa = state.H2OaElyte;
            
            params = model.elyteParams;
            con    = model.constants;

            F  = con.F;
            R  = con.R;
            c0 = params.c0;
            E0 = params.E0;

            state.Eelyte = E0 + R*T/(2*F)*log((H2Oa.^2).*(c0.^2).*(cOH.^-2));
            
        end

        function state = updateEinmr(model, state)

            T    = state.T;
            cOH  = state.cOHinmr;
            H2Oa = state.H2Oainmr;
            
            params = model.inmrParams;
            con    = model.constants;

            F  = con.F;
            R  = con.R;
            c0 = params.c0;
            E0 = params.E0;

            state.Einmr = E0 + R*T/(2*F)*log((H2Oa.^2).*(c0.^2).*(cOH.^-2));
            
        end

        function state = updateEtas(model, state)

            phiElyte = state.phiElyte;
            phiInmr  = state.phiInmr;
            Einmr    = state.Einmr;
            Eelyte   = state.Eelyte;

            state.etaElyte = E - phiElyte - Eelyte;
            state.etaInmr  = E - phiInmr - Einmr;
            
        end

        function state = updateReactionRateConstants(model, state)

            error('implemented in derived class')
            
        end
            
        function state = updateReactionRates(model, state)

            etaElyte = state.etaElyte;
            etaInmr  = state.etaInmr;
            T        = state.T;

            ir    = model.inmr;
            vsa   = model.volumetricSurfaceArea;
            Xinmr = model.Xinmr;
            alpha = model.alpha;
            n     = model.n;


            jElyte = (1 - X_inmr)*vsa*butlerVolmer(j0, alpha, n, etaElyte, T);
            jInmr = X_inmr*vsa*butlerVolmer(j0, alpha, n, etaInmr, T);

            state.elyteReactionRate = elyteReactionRate;
            state.inmrReactionRate  = inmrReactionRate;
            
        end


    end

end
