classdef LithiumPlatingLatz < BaseModel

    properties

        F = PhysicalConstants.F
        R = PhysicalConstants.R

        alphaPl    
        alphaStr   
        alphaChInt 

        kPl    
        kChInt 

        nPl0    
        muLiRef 

        nPlLimit 

        SEIFraction 
        MSEI        
        rhoSEI      
        deltaSEI0   
        sigmaSEI    

        useSEI 
        
    end

    methods

        function model = LithiumPlatingLatz(inputparams)
            model = model@BaseModel();
            fdnames = {'alphaPl'    , ...    
                       'alphaStr'   , ...   
                       'alphaChInt' , ... 
                       'kPl'        , ...    
                       'kChInt'     , ... 
                       'nPl0'       , ...    
                       'muLiRef'    , ... 
                       'nPlLimit'   , ... 
                       'SEIFraction', ... 
                       'MSEI'       , ...        
                       'rhoSEI'     , ...      
                       'deltaSEI0'  , ...   
                       'sigmaSEI'   , ...    
                       'useSEI'};
            model = dispatchParams(model, inputparams, fdnames);
            
        end

        function model = registerVarAndPropfuncNames(model)

            useSEI = model.useSEI;
            
            varnames = {};
            % potential of the solid electrode
            varnames{end + 1} = 'T';
            % potential of the solid electrode
            varnames{end + 1} = 'phiElectrode';
            %
            varnames{end + 1} = 'phiElectrolyte';
            %
            varnames{end + 1} = 'cElectrolyte';
            %
            varnames{end + 1} = 'nPl';
            %
            varnames{end + 1} = 'nPlAccum';
            %
            varnames{end + 1} = 'nPlCons';
            %
            varnames{end + 1} = 'platingFlux';
            %
            varnames{end + 1} = 'chemicalFlux';
            %
            varnames{end + 1} = 'etaPlating';
            %
            varnames{end + 1} = 'etaChemical';
            %
            varnames{end + 1} = 'activityPlated';
            %
            varnames{end + 1} = 'surfaceCoverage';

            if useSEI
                
                %
                varnames{end + 1} = 'nSEI';
                %
                varnames{end + 1} = 'nSEIAccum';
                %
                varnames{end + 1} = 'nSEICons';
                %
                varnames{end + 1} = 'SEIThickness';
                
            end
            
            model = model.registerVarNames(varnames);

            fn = @LithiumPlatingLatz.updateNPlAccum;
            fn = {fn, @(pf) PropFunction.accumFuncCallSetupFn(pf)};
            model = model.registerPropFunction({'nPlAccum', fn, {'nPl'}});

            fn = @LithiumPlatingLatz.updateActivityPlated;
            model = model.registerPropFunction({'activityPlated', fn, {'nPl'}});

            if useSEI
                fn = @LithiumPlatingLatz.updateEtaPlatingSEI;
                model = model.registerPropFunction({'etaPlating', fn, {'phiElectrode', 'phiElectrolyte', 'activityPlated', 'SEIThickness', 'T'}});
            else
                fn = @LithiumPlatingLatz.updateEtaPlating;
                model = model.registerPropFunction({'etaPlating', fn, {'phiElectrode', 'phiElectrolyte', 'activityPlated', 'T'}});
            end
            fn = @LithiumPlatingLatz.updateEtaChemical;
            model = model.registerPropFunction({'etaChemical', fn, {'activityPlated', 'T'}});

            fn = @LithiumPlatingLatz.updatePlatingFlux;
            model = model.registerPropFunction({'platingFlux', fn, {'cElectrolyte', 'etaPlating', 'T'}});

            fn = @LithiumPlatingLatz.updateChemicalFlux;
            model = model.registerPropFunction({'chemicalFlux', fn, {'etaChemical', 'T'}});

            fn = @LithiumPlatingLatz.updateNPlCons;
            model = model.registerPropFunction({'nPlCons', fn, {'nPlAccum', 'platingFlux', 'chemicalFlux'}});

            fn = @LithiumPlatingLatz.updateSurfaceCoverage;
            model = model.registerPropFunction({'surfaceCoverage', fn, {'nPl'}});

            if useSEI
                
                fn = @LithiumPlatingLatz.updateNSEIAccum;
                fn = {fn, @(pf) PropFunction.accumFuncCallSetupFn(pf)};
                model = model.registerPropFunction({'nSEIAccum', fn, {'nSEI'}});

                fn = @LithiumPlatingLatz.updateNSEICons;
                model = model.registerPropFunction({'nSEICons', fn, {'platingFlux', 'nSEIAccum'}});

                fn = @LithiumPlatingLatz.updateSEIThickness;
                model = model.registerPropFunction({'SEIThickness', fn, {'nSEI'}});

            end
        end

        function state = updateNPlAccum(model, state, state0, dt)
            state.nPlAccum = (state.nPl - state0.nPl) / dt;
        end

        function state = updateActivityPlated(model, state)
            
            nPl = state.nPl;
            n0 = model.nPl0;
            state.activityPlated = nPl ./ (nPl + n0);
            
        end

        function state = updateEtaPlatingSEI(model, state)
            
            phiS = state.phiElectrode;
            phiE = state.phiElectrolyte;
            aPl = state.activityPlated;

            RSEI = state.SEIThickness / model.sigmaSEI;
            j = model.kPl;  

            eta = phiS - phiE - j * model.F * RSEI + (model.R * state.T / model.F) .* log(aPl);
            state.etaPlating = eta;
        end

        function state = updateEtaPlating(model, state)
            
            phiS = state.phiElectrode;
            phiE = state.phiElectrolyte;
            aPl  = state.activityPlated;
            T    = state.T;

            eta = phiS - phiE  + (model.R * T / model.F) .* log(aPl);
            state.etaPlating = eta;
        end

        function state = updateEtaChemical(model, state)
            
            aPl = state.activityPlated;
            T   = state.T;
            
            eta = -(model.R * T / model.F) .* log(aPl);
            state.etaChemical = eta;
            
        end

        function state = updatePlatingFlux(model, state)

            R = model.R;
            F = model.F;
            
            eta = state.etaPlating;
            ce  = state.cElectrolyte;
            T   = state.T;
            
            i0 = model.kPl * ce.^model.alphaPl;
            j = i0 .* (exp((model.alphaPl * F * eta) / (R * T)) - ...
                       exp((-model.alphaStr * F * eta) / (R * T)));

            state.platingFlux = j ./ F;
            
        end

        function state = updateChemicalFlux(model, state)

            eta = state.etaChemical;
            T   = state.T;
            
            jCh = model.kChInt * (exp(0.5 * model.F * eta / (model.R * T)) - ...
                                  exp(-0.5 * model.F * eta / (model.R * T)));
            
            state.chemicalFlux = jCh ./ model.F;
            
        end

        function state = updateNPlCons(model, state)
            
            flux  = state.platingFlux - state.chemicalFlux;
            accum = state.nPlAccum;
            
            state.nPlCons = assembleConservationEquation(model, 0, 0, flux, accum);
            
        end

        function state = updateSurfaceCoverage(model, state)
            nPl = state.nPl;
            state.surfaceCoverage = min(nPl ./ model.nPlLimit, 1.0);
        end

        function state = updateNSEIAccum(model, state, state0, dt)
            state.nSEIAccum = (state.nSEI - state0.nSEI) / dt;
        end

        function state = updateNSEICons(model, state)
            
            flux  = state.platingFlux * model.SEIFraction;
            accum = state.nSEIAccum;
            
            state.nSEICons = assembleConservationEquation(model, 0, 0, flux, accum);
            
        end

        function state = updateSEIThickness(model, state)
            
            delta = model.deltaSEI0 + (model.MSEI * state.nSEI) / (model.rhoSEI * model.F);
            state.SEIThickness = delta;
        end

    end
end
