classdef ExchangeLayer < BaseModel

    properties

        inmrParams % struct with fields
        % kxch
        % cT
        % OH.z
        % kML
    end

    methods

        function model = ExchangeLayer(paramobj)

            model = model@BaseModel();
            
            fdnames = {'inmrParams'};
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

            %% Fluxes

            % flux in electrolyte
            varnames{end + 1} = 'H2OexchangeRate';
            % flux in ionomer
            varnames{end + 1} = 'OHexchangeRate';

            model = model.registerVarNames(varnames);

            % Assemble ionomer ion exchange rate
            fn = @() ExchangeLayer.updateOHexchange;
            inputnames = {'phiElyte', 'phiInmr', 'cOHelyte', 'cOHinmr'};
            model = model.registerPropFunction({'OHexchangeRate', fn, inputnames});

            % Assemble sorption rate
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

        %% TODO : fix this function
            eta      = state.eta;
            T        = state.T;
            cOHElyte = state.cOHElyte;

            ir = model.inmr;

            state.Rxch = ionomerExchange(ir.kxch, ir.cT, ir.OH.z, eta, T, cOHElyte);
        end

    end

end
