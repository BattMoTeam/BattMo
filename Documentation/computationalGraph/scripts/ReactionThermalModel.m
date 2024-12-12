classdef ReactionThermalModel < BaseModel

    properties

        Reaction
        Thermal
        
    end
    
    methods

        function model = ReactionThermalModel()

            model.Reaction = ReactionModel();
            model.Thermal  = ThermalModel();
            
        end
        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model

            model = registerVarAndPropfuncNames@BaseModel(model);

            fn = @ReactionThermalModel.updateOCP;
            inputnames = {{'Reaction', 'c_s'}, {'Thermal', 'T'}};
            model = model.registerPropFunction({{'Reaction', 'OCP'}, fn, inputnames});

            fn = @ReactionThermalModel.updateThermalSource;
            inputnames = {{'Reaction', 'R'}};
            model = model.registerPropFunction({{'Thermal', 'source'}, fn, inputnames});
            
        end

    end
    
end
