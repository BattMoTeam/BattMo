classdef EvolutionElectrode < BaseModel
    
    properties
        
        PorousTransportLayer
        CatalystLayer
        ExchangeLayer

        couplingTerms % Coupling terms
        couplingNames
        
    end

    methods

        function model = EvolutionElectrode(paramobj)

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

            %% Update coupling (electrode, ionomer) -> CatalystLayer

            fn = @() EvolutionElectrode.dispatchToCatalystAndExchangeLayers;

            phaseInd  = model.(ptl).phaseInd;
            liquidInd = model.(ptl).liquidInd;
            gasInd    = model.(ptl).gasInd;
            
            inputvarnames = {{ptl, 'phi'} , ...
                             VarName({ptl}, 'concentrations', liquidInd.ncomp, liquidInd.OH), ...
                             {ptl, 'H2Oa'}};
            outputvarnames = {'phiElyte', 'cOHelyte', 'H2OaElyte'};
            for iovar = 1 : numel(outputvarnames)
                outputvarname = {ctl, outputvarnames{iovar}};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                outputvarname = {exl, outputvarnames{iovar}};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            end

            inputvarnames = {VarName({ptl}, 'compGasPressures', gasInd.ncomp, gasInd.activeGas)};
            model = model.registerPropFunction({{ctl, 'pressureActiveGas'}, fn, inputvarnames});
            
            fn = @() EvolutionElectrode.updateSourceTerms;

            inputvarnames = {{ctl, 'elyteReactionRate'}, ...
                             {exl, 'OHexchangeRate'}};
            model = model.registerPropFunction({{ptl, 'OHSource'}, fn, inputvarnames});

            inputvarnames = {{ctl, 'elyteReactionRate'}, ...
                             {exl, 'H2OexchangeRate'}  , ...
                             {ptl, 'H2OliquidVaporExchangeRate'}};
            model = model.registerPropFunction({{ptl, 'H2OliquidSource'}, fn, inputvarnames});

            inputvarnames = {{ctl, 'elyteReactionRate'}, ...
                             {ctl, 'inmrReactionRate'}};
            gasInd = model.(ptl).gasInd;
            igas = gasInd.activeGas;
            ngas = gasInd.ncomp;
            model = model.registerPropFunction({VarName({ptl}, 'compGasSources', ngas, igas), fn, inputvarnames});


            fn = @() EvolutionElectrode.updateSourceTerms;

            inputvarnames = {{ctl, 'elyteReactionRate'}, ...
                             {exl, 'OHexchangeRate'}};
            model = model.registerPropFunction({{ptl, 'OHSource'}, fn, inputvarnames});
            
        end

        function state = dispatchTemperature(model, state)

            T = state.T;

            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            
            state.(ptl).T = T;
            state.(ctl).T = T;
            
        end
        
        function state = dispatchToCatalystAndExchangeLayers(model, state)

            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';
            
            lind = model.(ptl).liquidInd;
            gInd = model.(ptl).gasInd;

            coupterms = model.couplingTerms;
            coupnames = model.couplingNames;
            
            coupterm = getCoupTerm(coupterms, {ptl, ctl}, coupnames);
            coupcells = coupterm.couplingcells;
            
            phi  = state.(ptl).phi;
            cH2O = state.(ptl).concentrations{lind.H2Oliquid};
            cOH  = state.(ptl).concentrations{lind.OH};
            H2Oa = state.(ptl).H2Oactivity,
            pag  = state.(ptl).compGasPressures{gInd.activeGas};
            
            state.(ctl).phiElyte(coupcells(:, 2))          = phi(coupcells(:, 1));
            state.(ctl).cOHElyte(coupcells(:, 2))          = cOH(coupcells(:, 1));
            state.(ctl).cHElyte(coupcells(:, 2))           = cH(coupcells(:, 1));
            state.(ctl).H2OaElyte(coupcells(:, 2))         = H2Oa(coupcells(:, 1));
            state.(ctl).pressureActiveGas(coupcells(:, 2)) = pag(coupcells(:, 1));

            % coupling to ExchangeLayer
            coupterm = getCoupTerm(coupterms, {ptl, exl}, coupnames);
            coupcells = coupterm.couplingcells;
            state.(exl).phiElyte(coupcells(:, 2)) = phi(coupcells(:, 1));
            state.(exl).cOHElyte(coupcells(:, 2)) = cOH(coupcells(:, 1));
            state.(exl).H2aElyte(coupcells(:, 2)) = H2a(coupcells(:, 1));
            
        end

        function state = updateSourceTerms(model, state)
            
            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';
            
            leps = model.(inm).liquidVolumeFraction;
            gInd = model.(ptl).gasInd;
            
            coupterms = model.couplingTerms;
            coupnames = model.couplingNames;
            
            coupterm = getCoupTerm(coupterms, {ptl, ctl}, coupnames);
            coupcells = coupterm.couplingcells;
            
            jElyte    = state.(ctl).jElyte(coupcells(:, 2));
            j         = state.(ctl).j(coupcells(:, 2));
            Rxch      = state.(ctl).Rxch(coupcells(:, 2));
            Rsorption = state.(ctl).Rsorption(coupcells(:, 2));
            
            OHSource        = 2*jElyte./(n*F) + Rxch.*leps;
            H2OliquidSource = -Rsorption;
            compH2GasSource = -j/(n*F).*H2.MW;
            
            ccs = coupcells(:, 1);
            state.(ptl).OHSource(ccs)                = OHSource;
            state.(ptl).H2OliquidSource(ccs)         = H2OliquidSource;
            state.(ptl).compGasSources{gInd.H2}(ccs) = compH2GasSource;
            
        end
    
    
        function state = dispatchFromCatalystLayer(model, state)
            
            state = model.dispatchCatalystLayerToPorousTransportLayer(state);
            
        end


    end
    
    
end

