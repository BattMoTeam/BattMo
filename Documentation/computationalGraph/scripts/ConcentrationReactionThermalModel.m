classdef ConcentrationReactionThermalModel < BaseModel

    properties

        Masses
        Thermal
        
    end
    
    methods

        function model = ConcentrationReactionThermalModel()

            model.Masses  = ConcentrationReactionModel();
            model.Thermal = ThermalModel();
            
        end

        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model

            model = registerVarAndPropfuncNames@BaseModel(model);

            fn = @ConcentrationReactionThermalModel.updateOCP;
            inputnames = {{'Masses', 'Reaction', 'c_s'}, {'Thermal', 'T'}};
            model = model.registerPropFunction({{'Masses', 'Reaction', 'OCP'}, fn, inputnames});

            fn = @ConcentrationReactionThermalModel.updateThermalSource;
            inputnames = {{'Masses', 'Reaction', 'R'}};
            model = model.registerPropFunction({{'Thermal', 'source'}, fn, inputnames});
            
        end

    end
    
end
