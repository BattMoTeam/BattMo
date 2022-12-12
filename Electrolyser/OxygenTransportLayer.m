classdef OxygenPorousTransportLayer < PorousTransportLayer
    

    properties
        
    end
    
    methods
        
        function model = OxygenPorousTransportLayer(paramobj)

            model = model@PorousTransportLayer(paramobj);
            
            % add the  O2 component in the indexing structures            
            model.compInd.O2 = model.compInd.activeGas;
            model.gasInd.O2  = model.gasInd.activeGas;
        end

        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@PorousTransportLayer(model);

            phaseInd  = model.phaseInd;
            liquidInd = model.liquidInd;
            gasInd    = model.gasInd;
            compInd   = model.compInd;

            ncomp   = compInd.ncomp;
            ngas    = gasInd.ncomp;
            nph     = phaseInd.nphase;
            nliquid = liquidInd.ncomp;
            nmobph  = numel(phaseInd.mobile);

            model = model.registerVarName('O2rhoeps');

            % assemble gas pressure using ideal gas law
            fn = @() OxygenPorousTransportLayer.updateGasPressure;
            inputnames = {'O2rhoeps', 'H2Ogasrhoeps', 'T'};
            model = model.registerPropFunction({VarName({}, 'phasePressures', nph, phaseInd.gas), fn, inputnames});
            model = model.registerPropFunction({VarName({}, 'compGasPressures', ngas), fn, inputnames});
            
            fn = @() OxygenPorousTransportLayer.updateGasViscosity;
            inputnames = {'T'};
            model = model.registerPropFunction({VarName({}, 'viscosities', nph, phaseInd.gas), fn, inputnames});

            fn = @() OxygenPorousTransportLayer.updateO2Accum;
            functionCallSetupFn = @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction);
            fn = {fn, functionCallSetupFn};
            inputnames = {'O2rhoeps'};
            model = model.registerPropFunction({VarName({}, 'compGasAccums', ngas, gasInd.O2), fn, inputnames});
            
        end

        function state = updateGasPressure(model, state)

            gasInd   = model.gasInd;
            phaseInd = model.phaseInd;
            MWH2O    = model.sp.H2O.MW;
            MWO2     = model.sp.O2.MW;
            
            mO2  = state.O2rhoeps;
            mH2O = state.H2Ogasrhoeps;
            vf   = state.volumeFractions{model.phaseInd.gas};

            pO2  = R*T*(mO2./MWO2)./vf
            pH2O = R*T*(mH20./MWH2)./vf
            
            state.compGasPressure{gasInd.O2}   = pO2;
            state.compGasPressure{gasInd.H2O}  = pH2O
            state.phasePressures{phaseInd.gas} = pO2 + pH2O;
        end
        
        
        function state = updateGasViscosity(model, state)
            T = state.T;

            state.viscosities(model.phaseInd.gas) = (0.1971 + T * (0.0803 - 3.99e-5 * T)) * 1e-6;
        end
        
        
    end


end

