classdef ExchangeReaction < BaseModel

    properties

        constants % physical constants

        kxch % Exchange rate
        OH   % structure with field
             % - z : number of charge
        kML  % Ionomer sorption coefficient

    end

    methods

        function model = ExchangeReaction(paramobj)

            model = model@BaseModel();
            
            fdnames = {'kxch', ...
                       'OH'  , ...
                       'kML' };
            model = dispatchParams(model, paramobj, fdnames);

            model.constants = PhysicalConstants;
            
        end

        function model = registerVarAndPropfuncNames(model)

        %% Declaration of the Dynamical Variables and Function of the model
        % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {};
            % Electric Potential from electrolyte
            varnames{end + 1} = 'phiElyte';
            % Electric Potential from ionomer
            varnames{end + 1} = 'phiInmr';
            % concentration of OH from electrotyte
            varnames{end + 1} = 'cOHelyte';
            % concentration of OH from ionomer
            varnames{end + 1} = 'cOHinmr';
            % Water activity from electrolyte
            varnames{end + 1} = 'H2OaElyte';
            % Water activity from ionomer
            varnames{end + 1} = 'H2OaInmr';
            % Water activity from ionomer
            varnames{end + 1} = 'T';
            
            %% Fluxes

            % H2O exchange rate [mol*s^-1 m^-3]
            varnames{end + 1} = 'H2OexchangeRate';
            % OH exchange rate [mol*s^-1 m^-3]
            varnames{end + 1} = 'OHexchangeRate';

            model = model.registerVarNames(varnames);

            % Assemble ionomer ion exchange rate [mol*s^-1 m^-3]
            % (OH-)_inmr <->> (OH-)_elyte
            % Here, the direction of the reaction that is indicated by the repeated arrow symbol corresponds to a positive computed reaction rate
            fn = @() ExchangeReaction.updateOHexchange;

            inputnames = {'phiElyte', 'phiInmr', 'cOHelyte', 'cOHinmr', 'T'};
            model = model.registerPropFunction({'OHexchangeRate', fn, inputnames});

            % Assemble sorption rate [mol*s^-1 m^-3]
            % (H2O)_inmr <->> (H2O)_elyte
            % Here, the direction of the reaction that is indicated by the repeated arrow symbol corresponds to a positive computed reaction rate
            fn = @() ExchangeReaction.updateSorption;

            inputnames = {'H2OaElyte', 'H2OaInmr'};
            model = model.registerPropFunction({'H2OexchangeRate', fn, inputnames});


        end

        function state = updateSorption(model, state)

            kML = model.kML;
            con = model.constants;

            H2OaElyte = state.H2OaElyte;
            H2OaInmr  = state.H2OaInmr;
            T         = state.T;
            
            % (H2O)_inmr <->> (H2O)_elyte
            % Here, the direction of the reaction that is indicated by the repeated arrow symbol corresponds to a positive computed reaction rate
            state.H2OexchangeRate = kML*(con.R).*T.*(log(H2OaInmr) - log(H2OaElyte));
            
        end

        function state = updateOHexchange(model, state)
        %   Follows an approach from Stanislaw, Gerhardt, and Weber ECS Trans 2019,
        %   as well as Jiangjin Liu et al JES 2021.
        %   The equation for R approaches the equation for Donnan equilibrium if
        %   kxch is large or R approaches zero. So using a large kxch value should
        %   enforce Donnan equilibrium.

            
            kxch = model.kxch;
            z = model.OH.z;
            F = model.constants.F;
            R = model.constants.R;
            
            T        = state.T;
            cElyte   = state.cOHelyte;
            cInmr    = state.cOHinmr;
            phiElyte = state.phiElyte;
            phiInmr  = state.phiInmr;
            
            % (OH-)_inmr <->> (OH-)_elyte
            % Here, the direction of the reaction that is indicated by the repeated arrow symbol corresponds to a positive computed reaction rate
            state.OHexchangeRate = kxch.*(cInmr.*(exp(z*F*(phiInmr - phiElyte)./(R*T))) - cElyte);
            
        end

    end

end
