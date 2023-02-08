classdef EvolutionElectrode < BaseModel
    
    properties
        
        PorousTransportLayer
        CatalystLayer
        ExchangeLayer

        couplingTerm % Coupling term between ptl and ctl
        
    end

    methods

        function model = EvolutionElectrode(paramobj)

            fdnames = {'G',
                       'couplingTerm'};
            model = dispatchParams(model, paramobj, fdnames);            
            
            switch paramobj.porousTransportLayerType
              case 'Hydrogen'
                model.PorousTransportLayer = HydrogenPorousTransportLayer(paramobj.PorousTransportLayer);
              case 'Oxygen'
                model.PorousTransportLayer = OxygenPorousTransportLayer(paramobj.PorousTransportLayer);
              otherwise
                error('porousTransportLayerType not recognized')
            end
            
            switch paramobj.catalystLayerType
              case 'Iridium'
                model.CatalystLayer = IridiumCatalystLayer(paramobj.CatalystLayer);
              case 'Platinium'
                model.CatalystLayer = PlatiniumCatalystLayer(paramobj.CatalystLayer);
              otherwise
                error('catalystLayerType not recognized')
            end
            
            model.ExchangeLayer = ExchangeLayer(paramobj.ExchangeLayer);
            
        end

        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {'T'};
            model = model.registerVarNames(varnames);
            
            % shortcuts
            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';
        
            fn = @() EvolutionElectrode.dispatchTemperature;
            model = model.registerPropFunction({{ptl, 'T'}, fn, {'T'}});
            model = model.registerPropFunction({{ctl, 'T'}, fn, {'T'}});
            model = model.registerPropFunction({{exl, 'T'}, fn, {'T'}});
            
            %% Update coupling (electrode, ionomer) -> CatalystLayer

            fn = @() EvolutionElectrode.dispatchToCatalystAndExchangeLayers;

            phaseInd  = model.(ptl).phaseInd;
            liquidInd = model.(ptl).liquidInd;
            gasInd    = model.(ptl).gasInd;
            
            inputvarnames = {{ptl, 'phi'}                                                     , ...
                             VarName({ptl}, 'concentrations', liquidInd.nliquid, liquidInd.OH), ...
                             VarName({ptl}, 'compGasPressures', gasInd.ngas, gasInd.activeGas), ...
                             {ptl, 'H2Oa'}                                                    , ...
                             {ptl, 'E'}};
            outputvarnames = {'phiElyte', 'cOHelyte', 'H2OaElyte'};
            for iovar = 1 : numel(outputvarnames)
                outputvarname = {ctl, outputvarnames{iovar}};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                outputvarname = {exl, outputvarnames{iovar}};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            end
            model = model.registerPropFunction({{ctl, 'E'}, fn, inputvarnames});
            model = model.registerPropFunction({{ctl, 'pressureActiveGas'}, fn, inputvarnames});
            
            fn = @() EvolutionElectrode.updateSourceTerms;

            inputvarnames = {{ctl, 'elyteOHsource'}             , ...
                             {exl, 'OHexchangeRate'}            , ...
                             {ctl, 'elyteH2Osource'}            , ...
                             {exl, 'H2OexchangeRate'}           , ...
                             {ptl, 'H2OvaporLiquidExchangeRate'}, ...
                             {ctl, 'activeGasSource'}};
            model = model.registerPropFunction({{ptl, 'OHsource'}, fn, inputvarnames});
            model = model.registerPropFunction({{ptl, 'H2OliquidSource'}, fn, inputvarnames});
            gasInd = model.(ptl).gasInd;
            igas = gasInd.activeGas;
            ngas = gasInd.ngas;
            model = model.registerPropFunction({VarName({ptl}, 'compGasSources', ngas, igas), fn, inputvarnames});
            
        end

        function state = dispatchTemperature(model, state)

            T = state.T;

            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';
            
            coupterm = model.couplingTerm;
            coupcells = coupterm.couplingcells;
            
            state.(ptl).T = T;
            state.(ctl).T(coupcells(:, 2), 1) = T(coupcells(:, 1));
            state.(exl).T(coupcells(:, 2), 1) = T(coupcells(:, 1));
            
        end
        
        function state = dispatchToCatalystAndExchangeLayers(model, state)

            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';
            
            lind = model.(ptl).liquidInd;
            gInd = model.(ptl).gasInd;

            coupterm = model.couplingTerm;
            coupcells = coupterm.couplingcells;
            
            phi  = state.(ptl).phi;
            cH2O = state.(ptl).concentrations{lind.H2O};
            cOH  = state.(ptl).concentrations{lind.OH};
            H2Oa = state.(ptl).H2Oa;
            pag  = state.(ptl).compGasPressures{gInd.activeGas};
            E    = state.(ptl).E;

            % initialization of the variables (takes care of AD)
            nc = model.(ctl).G.cells.num;
            initval = nan(nc, 1);
            [adsample, isAD] = getSampleAD(phi);
            if isAD
                initval = adsample.convertDouble(initval);
            end
                
            varnames = {'phiElyte', 'cOHelyte', 'H2OaElyte'};
            layers = {ctl, exl};
            for ivar = 1 : numel(varnames)
                varname = varnames{ivar};
                for ilayer = 1 : numel(layers)
                    layer = layers{ilayer};
                    state.(layer).(varname) = initval;
                end
            end
            state.(ctl).pressureActiveGas = initval;
            
            state.(ctl).phiElyte(coupcells(:, 2))          = phi(coupcells(:, 1));
            state.(ctl).cOHelyte(coupcells(:, 2))          = cOH(coupcells(:, 1));
            % state.(ctl).cHElyte(coupcells(:, 2))         = cH(coupcells(:, 1));
            state.(ctl).H2OaElyte(coupcells(:, 2))         = H2Oa(coupcells(:, 1));
            state.(ctl).pressureActiveGas(coupcells(:, 2)) = pag(coupcells(:, 1));
            
            state.(ctl).E = E;
            
            % coupling to ExchangeLayer
            state.(exl).phiElyte(coupcells(:, 2))  = phi(coupcells(:, 1));
            state.(exl).cOHelyte(coupcells(:, 2))  = cOH(coupcells(:, 1));
            state.(exl).H2OaElyte(coupcells(:, 2)) = H2Oa(coupcells(:, 1));
            
        end

        function state = updateSourceTerms(model, state)
            
            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';
            
            gInd     = model.(ptl).gasInd;
            coupterm = model.couplingTerm;
            gasMW    = model.(ptl).gasMW;
            vols     = model.(ptl).G.cells.volumes;
            
            coupcells = coupterm.couplingcells;

            H2OvlR = state.(ptl).H2OvaporLiquidExchangeRate;

            H2OexchR                  = 0*H2OvlR; % initialization so that get we get right AD and dimension
            H2OexchR(coupcells(:, 1)) = state.(exl).H2OexchangeRate(coupcells(:, 2));
            OHexchR                   = 0*H2OvlR; % initialization so that get we get right AD and dimension
            OHexchR(coupcells(:, 1))  = state.(exl).OHexchangeRate(coupcells(:, 2));
            
            elyteOHsource                    = 0*H2OvlR; % initialization so that get we get right AD and dimension
            elyteOHsource(coupcells(:, 1))   = state.(ctl).elyteOHsource(coupcells(:, 2));
            activeGasSource                  = 0*H2OvlR; % initialization so that get we get right AD and dimension
            activeGasSource(coupcells(:, 1)) = state.(ctl).activeGasSource(coupcells(:, 2));
            elyteH2Osource                   = 0*H2OvlR; % initialization so that get we get right AD and dimension
            elyteH2Osource(coupcells(:, 1))  = state.(ctl).elyteH2Osource(coupcells(:, 2));
            
            state.(ptl).OHsource                       = vols.*(elyteOHsource + OHexchR);
            state.(ptl).H2OliquidSource                = vols.*(elyteH2Osource + H2OexchR + H2OvlR);
            state.(ptl).compGasSources{gInd.activeGas} = gasMW.*vols.*activeGasSource;

        end
    
    
        function state = dispatchFromCatalystLayer(model, state)
            
            state = model.dispatchCatalystLayerToPorousTransportLayer(state);
            
        end


    end
    
    
end

