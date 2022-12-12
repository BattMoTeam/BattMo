classdef EvolutionElectrode < BaseModel
    
    properties
        
        PorousTransportLayer
        CatalystLayer
        ExchangeLayer
        
    end

    methods

        function model = EvolutionElectrode(paramobj)

            model.PorousTransportLayer = PorousTransportLayer(paramobj.PorousTransportLayer);
            model.CatalystLayer        = CatalystLayer(paramobj.CatalystLayer);
            model.ExchangeLayer        = ExchangeLayer(paramobj.ExchangeLayer);
            
        end

        
        function model = registerVarAndPropfuncNames(model)

            varnames = {'T'};
            model = registerVarAndPropfuncNames@BaseModel(model);

            % shortcuts
            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';
        
            fn = @() Electrolyser.dispatchTemperature;
            model = model.registerPropFunction({{ptl, 'T'}, fn, {'T'}});
            model = model.registerPropFunction({{ctl, 'T'}, fn, {'T'}});

            %% Update coupling (electrode, ionomer) -> CatalystLayer

            fn = @() Electrolyser.dispatchToCatalystLayersAndInterfaces;

            phaseInd  = model.(ptl).phaseInd;
            liquidInd = model.(ptl).liquidInd;
            inputvarnames = {{ptl, 'phi'} , ...
                             VarName({ptl}, 'concentrations', liquidInd.ncomp, liquidInd.OH), ...
                             {ptl, 'H2Oa'}};
            outputvarnames = {'phiElyte', 'cOHelyte', 'H2OaElyte'};
            for iovar = 1 : numel(outputvarnames)
                outputvarname = {cat, outputvarnames{iovar}};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                outputvarname = {exl, outputvarnames{iovar}};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            end

            inputvarnames = {VarName({elde}, 'phasePressures', phaseInd.nphase, phaseInd.gas)};
            model = model.registerPropFunction({{cat, 'pressureActiveGas'}, fn, inputvarnames});
            
            fn = @() Electrolyser.updateSourceTerms;

            inputvarnames = {{cat, 'elyteReactionRate'}, ...
                             {exl, 'OHexchangeRate'}};
            model = model.registerPropFunction({{ptl, 'OHSource'}, fn, inputvarnames});

            inputvarnames = {{cat, 'elyteReactionRate'}, ...
                             {exl, 'H2OexchangeRate'}  , ...
                             {exl, 'H2OliquidVaporExchangeRate'}};
            model = model.registerPropFunction({{ptl, 'H2OliquidSource'}, fn, inputvarnames});

            inputvarnames = {{cat, 'elyteReactionRate'}, ...
                             {cat, 'inmrReactionRate'}};
            gasInd = model.(elde).gasInd;
            igas = gasInd.activeGas;
            ngas = gasInd.ncomp;
            model = model.registerPropFunction({VarName({ptl}, 'compGasSources', ngas, igas), fn, inputvarnames});

        end

        function state = dispatchTemperature(model, state)

            T = state.T;

            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            
            state.(ptl).T = T;
            state.(ctl).T = T;
            
        end
        
        function state = dispatchToCatalystLayersAndInterfaces(model, state)

            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';
            
            coupterms = model.couplingTerms;
            coupnames = model.couplingNames;
            
            coupterm = getCoupTerm(coupterms, {ptl, ctl}, coupnames);
            coupcells = coupterm.couplingcells;
            state.(ctl).phiElyte(coupcells(:, 2)) = state.(ptl).phi(coupcells(:, 1));
            
            gind = model.(ptl).phaseInd.gas;
            m  = state.(ptl).compGasMasses;
            T  = state.(ptl).T;
            vg = state.(ptl).volumeFraction{gind};
            
            gInd = model.(ptl).gasInd;
            
            % compute mol amount for each gas component

            mH2O = m{gInd.H2Ogas};
            mActGas = m{gInd.activeGas}
            
            mH2O    = mH2O(coupcells(:, 1));
            mActGas = mActGas(coupcells(:, 1));
            vg      = vg(coupcells(:, 1));
            T       = T(coupcells(:, 1));
            
            p{gInd.H2Ogas} = mH2O./(model.(ptl).H2O.MW).*R.*T./vg;
            p{gInd.activeGas} = mActGas./(model.(ptl).activeGas.MW).*R.*T./vg;
            
            state.(ctl).pressureActiveGas(coupcells(:, 2)) = p;





            
            coupterm = getCoupTerm(coupterms, {ptl, exl}, coupnames);
            coupcells = coupterm.couplingcells;
            state.(exl).phiElyte(coupcells(:, 2)) = state.(ptl).phi(coupcells(:, 1));
            
            
        end
        
        function state = dispatchFromCatalystLayer(model, state)
            
            state = model.dispatchCatalystLayerToPorousTransportLayer(state);
            
        end


    end
    
    
end

