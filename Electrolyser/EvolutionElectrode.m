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
                MWs(1) = paramobj.PorousTransportLayer.sp.H2O.MW;
                MWs(2) = paramobj.PorousTransportLayer.sp.H2.MW;
                paramobj.PorousTransportLayer.Boundary.MWs = MWs;
                model.PorousTransportLayer = HydrogenPorousTransportLayer(paramobj.PorousTransportLayer);
              case 'Oxygen'
                MWs(1) = paramobj.PorousTransportLayer.sp.H2O.MW;
                MWs(2) = paramobj.PorousTransportLayer.sp.O2.MW;
                paramobj.PorousTransportLayer.Boundary.MWs = MWs;
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
            
            inputvarnames = {{ptl, 'phi'}                                                      , ...
                             VarName({ptl}, 'concentrations', liquidInd.ncomp, liquidInd.OH)   , ...
                             VarName({ptl}, 'compGasPressures', gasInd.ncomp, gasInd.activeGas), ...
                             {ptl, 'H2Oa'}};
            outputvarnames = {'phiElyte', 'cOHelyte', 'H2OaElyte'};
            for iovar = 1 : numel(outputvarnames)
                outputvarname = {ctl, outputvarnames{iovar}};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                outputvarname = {exl, outputvarnames{iovar}};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            end
            model = model.registerPropFunction({{ctl, 'pressureActiveGas'}, fn, inputvarnames});
            
            fn = @() EvolutionElectrode.updateSourceTerms;

            inputvarnames = {{ctl, 'elyteOHsource'}, ...
                             {exl, 'OHexchangeRate'}};
            model = model.registerPropFunction({{ptl, 'OHSource'}, fn, inputvarnames});

            inputvarnames = {{ctl, 'elyteH2Osource'}, ...
                             {exl, 'H2OexchangeRate'}  , ...
                             {ptl, 'H2OvaporLiquidExchangeRate'}};
            model = model.registerPropFunction({{ptl, 'H2OliquidSource'}, fn, inputvarnames});

            inputvarnames = {{ctl, 'activeGasSource'}};
            gasInd = model.(ptl).gasInd;
            igas = gasInd.activeGas;
            ngas = gasInd.ncomp;
            model = model.registerPropFunction({VarName({ptl}, 'compGasSources', ngas, igas), fn, inputvarnames});
            
        end

        function state = dispatchTemperature(model, state)

            T = state.T;

            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';
            
            state.(ptl).T = T;
            state.(ctl).T = T;
            state.(exl).T = T;
            
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
            cH2O = state.(ptl).concentrations{lind.H2Oliquid};
            cOH  = state.(ptl).concentrations{lind.OH};
            H2Oa = state.(ptl).H2Oa;
            pag  = state.(ptl).compGasPressures{gInd.activeGas};
            
            state.(ctl).phiElyte(coupcells(:, 2))          = phi(coupcells(:, 1));
            state.(ctl).cOHElyte(coupcells(:, 2))          = cOH(coupcells(:, 1));
            % state.(ctl).cHElyte(coupcells(:, 2))           = cH(coupcells(:, 1));
            state.(ctl).H2OaElyte(coupcells(:, 2))         = H2Oa(coupcells(:, 1));
            state.(ctl).pressureActiveGas(coupcells(:, 2)) = pag(coupcells(:, 1));

            % coupling to ExchangeLayer
            state.(exl).phiElyte(coupcells(:, 2))  = phi(coupcells(:, 1));
            state.(exl).cOHElyte(coupcells(:, 2))  = cOH(coupcells(:, 1));
            state.(exl).H2OaElyte(coupcells(:, 2)) = H2Oa(coupcells(:, 1));
            
        end

        function state = updateSourceTerms(model, state)
            
            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';
            
            lvf       = model.(inm).liquidVolumeFraction;
            gInd      = model.(ptl).gasInd;
            coupterms = model.couplingTerms;
            coupnames = model.couplingNames;
            gasMW     = model.(ptl).gasMW;
            vols      = model.(ptl).G.cells.volumes;
            
            coupterm = getCoupTerm(coupterms, {ptl, ctl}, coupnames);
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
            
            state.(ptl).OHSource                       = vols.*(elyteOHsource + OHexchR);
            state.(ptl).H2OliquidSource                = vols.*(elyteH2Osource + H2OexchR + H2OvlR);
            state.(ptl).compGasSources{gInd.activeGas} = gasMW.*vols.*activeGasSource;

        end
    
    
        function state = dispatchFromCatalystLayer(model, state)
            
            state = model.dispatchCatalystLayerToPorousTransportLayer(state);
            
        end


    end
    
    
end

