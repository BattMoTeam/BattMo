classdef ExchangeLayer < BaseModel

    properties

        kxch %
        OH % structure with field
           % - z : number of charge
        kML % 
        
    end

    methods

        function model = ExchangeLayer(paramobj)

            model = model@BaseModel();
            
            fdnames = {'kxch', ...
                       'IH'  , ...
                       'KML' };
            model = dispatchParams(model, paramobj, fdnames);
            
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
            % (OH-)_inmr <-> (OH-)_elyte
            fn = @() ExchangeLayer.updateOHexchange;
            inputnames = {'phiElyte', 'phiInmr', 'cOHelyte', 'cOHinmr', 'T'};
            model = model.registerPropFunction({'OHexchangeRate', fn, inputnames});

            % Assemble sorption rate [mol*s^-1 m^-3]
            % (H2O)_inmr <-> (H2O)_elyte  
            fn = @() ExchangeLayer.updateSorption;
            inputnames = {'H2OaElyte', 'H2OaInmr'};
            model = model.registerPropFunction({'H2OexchangeRate', fn, inputnames});


        end

        function state = updateSorption(model, state)

            H2OaElyte = state.H2OaElyte;
            H2OaInmr = state.H2OaInmr;

            kML = model.inmr.kML;

            state.Rsorption = ionomerSorption(kML, H2OaElyte, H2OaInmr);
        end

        function state = updateOHexchange(model, state)

            kxch = model.kxch;
            
            eta      = state.eta;
            T        = state.T;
            cOHElyte = state.cOHElyte;
            cOHInmr  = state.cOHInmr;
            
            state.Rxch = ionomerExchange(kxch, cOHinmr, OH.z, eta, T, cOHElyte);
            
        end

    end

end
