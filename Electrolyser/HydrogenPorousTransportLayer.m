classdef HydrogenPorousTransportLayer < PorousTransportLayer

    properties
        
    end
    
    methods
        
        function model = HydrogenPorousTransportLayer(paramobj)

            model = model@PorousTransportLayer(paramobj);

            % add the H2 component in the indexing structures
            model.compInd.H2 = model.compInd.activeGas;
            model.gasInd.H2  = model.gasInd.activeGas;
            
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

            model = model.registerVarName('H2rhoeps');

            % assemble gas pressure using ideal gas law
            fn = @() HydrogenPorousTransportLayer.updateGasPressure;
            inputnames = {'H2rhoeps', 'H2Ogasrhoeps', 'T'};
            model = model.registerPropFunction({VarName({}, 'phasePressures', nph, phaseInd.gas), fn, inputnames});
            model = model.registerPropFunction({VarName({}, 'compGasPressures', ngas), fn, inputnames});

            % assemble H2mass
            fn = @() HydrogenPorousTransportLayer.updateH2GasMass;
            inputnames = {'H2rhoeps'};
            model = model.registerPropFunction({VarName({}, 'compGasMasses', ngas, gasInd.H2), fn, inputnames});
            
            fn = @() HydrogenPorousTransportLayer.updateGasViscosity;
            inputnames = {'T'};
            model = model.registerPropFunction({VarName({}, 'viscosities', nph, phaseInd.gas), fn, inputnames});
            
            
        end

        function state = updateH2GasMass(model, state)

            gasInd = model.gasInd;

            state.compGassMasses{gasInd.H2} = state.H2rhoeps;

        end
        
        function state = updateGasPressure(model, state)

            gasInd   = model.gasInd;
            phaseInd = model.phaseInd;
            MWH2O    = model.sp.H2O.MW;
            MWH2     = model.sp.H2.MW;
            
            mH2  = state.H2rhoeps;
            mH2O = state.H2Ogasrhoeps;
            vf   = state.volumeFractions{model.phaseInd.gas};

            pH2  = R*T*(mH2./MWH2)./vf
            pH2O = R*T*(mH20./MWH2)./vf
            
            state.compGasPressure{gasInd.H2}   = pH2;
            state.compGasPressure{gasInd.H2O}  = pH2O
            state.phasePressures{phaseInd.gas} = pH2 + pH2O;
            
        end
        
        
        function state = updateGasViscosity(model, state)
            
            T = state.T;
            warning('check that one! this was taken directly from OxygenPorousTransportLayer');
            state.viscosities(model.mobPhaseInd.gas) = (0.1971 + T * (0.0803 - 3.99e-5 * T)) * 1e-6;
            
        end
        
        
    end


end


