classdef OxygenElectrode < OpenElectrode
    

    properties
        
    end
    
    methods
        
        function model = OxygenElectrode(paramobj)

            model = model@OpenElectrode(paramobj);
            
            % add the  O2 component in the indexing structures            
            model.compInd.O2  = 5;
            model.phaseInd.O2 = 2;
            model.gasInd.O2   = 2;            
        end

        
        function model = registerVarAndPropfuncNames(model)
            
        %% Declaration of the Dynamical Variables and Function of the model
        % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@OpenElectrode(model);

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
            fn = @() HydrogenElectrode.updateGasPressure;
            inputnames = {'O2rhoeps', 'H2Ogasrhoeps', 'T'};
            model = model.registerPropFunction({VarName({}, 'pressures', nph, phaseInd.gas), fn, inputnames});
            
            if model.useZeroDmodel
                % assemble mass of O2
                fn = @() HydrogenElectrode.updatemassO20;
                inputnames = {VarName({}, 'compGasBcSources', ngas, gasInd.O2), ...
                              VarName({}, 'compGasSources', ngas, gasInd.O2), ...
                              VarName({}, 'accumTerms', ncomp, compInd.O2)};
                model = model.registerPropFunction({VarName({}, 'gasMassCons', nph, gasInd.O2), fn, inputnames});

                fn = @() OpenElectrode.updateO2Accum;
                inputnames = {'O2rhoeps'};
                functionCallSetupFn = @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction);
                model = model.registerPropFunction({VarName({}, 'accumTerms', ncomp, compInd.O2), fn, inputnames, functionCallSetupFn});
                
            else

                error('not yet implemented');
                
            end
            
        end

        function state = updateGasPressure(model, state)
            
            MWH2O = model.sp.H2O.MW;
            MWO2 = model.sp.O2.MW;
            
            mO2  = state.O2rhoeps;
            mH2O = state.H2Ogasrhoeps;
            vfg  = state.volumeFractions{model.phaseInd.gas};
            

            state.pressures{model.phaseInd.gas} = (mH20./MWH2O + mO2./MWO2).*R.*T./vfg;
        end
        
        
        function state = updateGasViscosity(model, state)
            T = state.T;

            state.viscosities(model.phaseInd.gas) = (0.1971 + T * (0.0803 - 3.99e-5 * T)) * 1e-6;
        end
        
        
    end


end

