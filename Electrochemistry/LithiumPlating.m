classdef LithiumPlating < BaseModel

    properties

        alphaC = 0.5
        alphaA = 0.5
        kLi = 1e-9
        
    end

    methods

        function model = LithiumPlating(inputparams)

            model = model@BaseModel();
            
            % fdnames = {'alphaC', ...
            %            'alphaA', ...
            %            'kLi'};

            % model = dispatchParams(model, inputparams, fdnames);
            
        end

        function model = registerVarAndPropfuncNames(model)

            varnames = {'phiElectrode'  , ...
                        'phiElectrolyte', ...
                        'cElectrolyte'  , ...
                        'concentration' , ...
                        'flux'          , ...
                        'massAccum'     , ...
                        'massCons'};

            model = model.registerVarNames(varnames);
            
            fn = @LithiumPlating.updateMassCons; 
            model = model.registerPropFunction({'massCons', fn, {'massAccum', 'flux'}});
            
            fn = @LithiumPlating.updateMassAccum; 
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({'massAccum', fn, {'concentration'}});

            fn = @LithiumPlating.updateFlux; 
            model = model.registerPropFunction({'flux', fn, {'phiElectrode', 'phiElectrolyte', 'cElectrolyte', 'concentration'}});
            
        end
    end
    
end

