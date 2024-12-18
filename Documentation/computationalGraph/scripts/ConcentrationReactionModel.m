classdef ConcentrationReactionModel < BaseModel

    properties

        Reaction
        Solid
        Elyte
        
    end
    
    methods

        function model = ConcentrationReactionModel()

            model.Reaction = ReactionModel();
            model.Solid    = ConcentrationModel();
            model.Elyte    = ConcentrationModel();
            
        end
        
        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model

            model = registerVarAndPropfuncNames@BaseModel(model);

            fn = @ConcentrationReactionModel.updateConcentrationSource;
            inputnames = {{'Reaction', 'R'}};
            model = model.registerPropFunction({{'Solid', 'source'}, fn, inputnames});
            model = model.registerPropFunction({{'Elyte', 'source'}, fn, inputnames});
            
            fn = @ConcentrationReactionModel.updateReactionConcentrationE;
            inputnames = {{'Elyte', 'c'}};
            model = model.registerPropFunction({{'Reaction', 'c_e'}, fn, inputnames});

            fn = @ConcentrationReactionModel.updateReactionConcentrationS;
            inputnames = {{'Solid', 'c'}};
            model = model.registerPropFunction({{'Reaction', 'c_s'}, fn, inputnames});

        end

    end
    
end
