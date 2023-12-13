classdef SeaWaterElectrode < ElectronicComponent

    properties
        V                    % Molar Volume
        externalCouplingTerm % 
    end
    
    methods

        function model = SeaWaterElectrode(inputparams)
            
            model = model@ElectronicComponent(inputparams);
            
            fdnames = {'V', ...
                       'externalCouplingTerm'};
            model = dispatchParams(model, inputparams, fdnames);
            
            model.operators = localSetupOperators(model.G);
            
        end
        
        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames@ElectronicComponent(model);
            
            varnames = {};
            % Volume fraction
            varnames{end + 1} = 'volumeFraction';
            % Accumulation term for active component
            varnames{end + 1} = 'accumTerm';
            % Source term for mass conservation equation
            varnames{end + 1} = 'sourceTerm';
            % Mass conservation equation
            varnames{end + 1} = 'massCons';
            
            model = model.registerVarNames(varnames);
            
            fn = @() SeaWaterElectrode.updateMassCons;
            inputnames = {'accumTerm', 'sourceTerm'};
            model = model.registerPropFunction({'massCons', fn, inputnames});
            
            fn = @() SeaWaterElectrode.updateConductivity;
            inputnames = {'volumeFraction'};
            model = model.registerPropFunction({'conductivity', fn, inputnames});
            
            fn = @() SeaWaterElectrode.updateAccumTerm;
            fn = {fn, @(propfunc) PropFunction.accumFuncCallSetupFn(propfunc)};
            inputnames = {'volumeFraction'};
            model = model.registerPropFunction({'accumTerm', fn, inputnames});
            
        end

        function state = updateAccumTerm(model, state, state0, dt)

            vols = model.G.cells.volumes;
            state.accumTerm = 1/dt*vols.*(state.volumeFraction - state0.volumeFraction);
            
        end


    
        function state = updateConductivity(model, state)

            sigma = model.EffectiveElectricalConductivity; % this prefix "Effective" should be changed in this property name
            
            vf = state.volumeFraction;
            
            state.conductivity = (vf.^1.5).*sigma;
            
        end
        
        
        function state = updateMassCons(model, state)
            accumTerm = state.accumTerm;
            sourceTerm = state.sourceTerm;
            
            state.massCons = accumTerm - sourceTerm;
        end
        
    end
    
end
