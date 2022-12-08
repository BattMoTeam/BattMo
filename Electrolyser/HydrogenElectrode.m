classdef HydrogenElectrode < OpenElectrode

    properties
        
    end
    
    methods
        
        function model = HydrogenElectrode(paramobj)

            model = model@OpenElectrode(paramobj);

            % add the H2 component in the indexing structures
            model.compInd.H2  = 5;
            model.phaseInd.H2 = 2;
            model.gasInd.H2   = 2;
            
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

            model = model.registerVarName('H2rhoeps');

            % assemble gas pressure using ideal gas law
            fn = @() HydrogenElectrode.updateGasPressure;
            inputnames = {'H2rhoeps', 'H2Ogasrhoeps', 'T'};
            model = model.registerPropFunction({VarName({}, 'pressures', nph, phaseInd.gas), fn, inputnames});

            fn = @() HydrogenElectrode.updateGasViscosity;
            inputnames = {'T'};
            model = model.registerPropFunction({VarName({}, 'viscosities', nph, phaseInd.gas), fn, inputnames});
            
            fn = @() HydrogenElectrode.updateH2Accum;
            functionCallSetupFn = @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction);
            fn = {fn, functionCallSetupFn};
            inputnames = {'H2rhoeps'};
            model = model.registerPropFunction({VarName({}, 'compGasAccums', ngas, gasInd.H2), fn, inputnames});
            
        end
        
        
    end


end


