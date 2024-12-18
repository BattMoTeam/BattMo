classdef ConcentrationModel < BaseModel

    methods

        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            % Concentration
            varnames{end + 1} = 'c';
            % Mass accumulation term
            varnames{end + 1} = 'massAccum';
            % Source term
            varnames{end + 1} = 'source';
            % Mass conservation equation
            varnames{end + 1} = 'massCons';
            
            model = model.registerVarNames(varnames);
            
            fn = @ConcentrationModel.updateMassAccum;
            inputnames = {'c'};
            model = model.registerPropFunction({'massAccum', fn, inputnames});

            fn = @ConcentrationModel.updateMassCons;
            inputnames = {'massAccum', 'source'};
            model = model.registerPropFunction({'massCons', fn, inputnames});
            
        end

    end
    
end
