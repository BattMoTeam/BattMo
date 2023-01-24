classdef CatalystLayer < BaseModel

    properties

        constants

        compInd   % mapping structure for component indices
        phaseInd  % mapping structure for phase indices
        liquidInd % mapping structure for component indices
        gasInd    % mapping structure for component indices

        j0 % Exchange current density
        E0
        sp % species struct with field
        % - OH.z  : Charge
        % - OH.c0 : OH reference concentration
        
        Xinmr % Fraction of specific area that is coversed with ionomer [-]
        volumetricSurfaceArea % Volumetric surface area [m^ -1]
        
    end

    methods

        function model = CatalystLayer(paramobj)
            
            model = model@BaseModel();
            
            fdnames = { 'j0'   , ...
                        'E0'   , ...
                        'sp'   , ...
                        'Xinmr', ...
                        'volumetricSurfaceArea'};
            
            model = dispatchParams(model, paramobj, fdnames);

            model.constants = PhysicalConstants();
            
        end

        function model = registerVarAndPropfuncNames(model)

        %% Declaration of the Dynamical Variables and Function of the model
        % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            names = {};
            % Temperature
            varnames = {'T'};
            % CatalystLayer electrical potential (single value)
            varnames{end + 1} = {'E'};
            % CatalystLayer total current (scalar value)
            varnames{end + 1} = {'I'};
            % Electric Potential in electrolyte and ionomer
            varnames{end + 1} = 'phiElyte';
            varnames{end + 1} = 'phiInmr';
            % concentration of OH in electrotyte and ionomer
            varnames{end + 1} = 'cOHelyte';
            varnames{end + 1} = 'cOHinmr';
            % Equilibrium Potential for electrolyte and ionomer
            varnames{end + 1} = 'Eelyte';
            varnames{end + 1} = 'Einmr';
            % Partial pressure of the active gas (H2 or O2 for exampler) in the porous transport layer
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
            
            % current source [A]
            varnames{end + 1} = 'eSource'; 
            % The following source terms are per volume. The units are [mol s^-1 m^-3]
            varnames{end + 1} = 'activeGasSource';
            varnames{end + 1} = 'elyteH2Osource';
            varnames{end + 1} = 'elyteOHsource';
            varnames{end + 1} = 'inmrOHsource';
            
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
            inputnames = {'phiElyte', 'E' 'Eelyte'};
            model = model.registerPropFunction({'etaElyte', fn, inputnames});
            inputnames = {'phiInmr', 'E' 'Einmr'};
            model = model.registerPropFunction({'etaInmr', fn, inputnames});            

            % update source terms
            fn = @() CatalystLayer.updateSources;
            inputnames = {'inmrReactionRate'};
            model = model.registerPropFunction({'inmrOHsource', fn, inputnames});
            inputnames = {'elyteReactionRate'};
            model = model.registerPropFunction({'elyteOHsource', fn, inputnames});
            inputnames = {'elyteReactionRate', 'inmrReactionRate'};
            model = model.registerPropFunction({'elyteH2Osource', fn, inputnames});
            model = model.registerPropFunction({'activeGasSource', fn, inputnames});
            model = model.registerPropFunction({'eSource', fn, inputnames});
            
            % Assemble the reaction rates
            fn = @() CatalystLayer.updateReactionRates;
            inputnames = {'elyteReactionRateConstant', 'etaElyte'};
            model = model.registerPropFunction({'elyteReactionRate', fn, inputnames});
            inputnames = {'inmrReactionRateConstant', 'etaInmr'};
            model = model.registerPropFunction({'inmrReactionRate', fn, inputnames});            

            fn = @() CatalystLayer.updateI;
            inputnames = {'eSource'};
            model = model.registerPropFunction({'I', fn, inputnames});

        end

        function state = updateSources(model, state)
            error('virtual function. Should be overloaded by specific catalystlayer');
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

        function state = updateEtas(model, state)

            E        = state.E;
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

        function state = updateI(model, state)

            state.I = sum(state.eSource);
            
        end
        

    end

end
