classdef UncoupledReactionThermalModel < BaseModel
    properties
        Reaction
        Thermal
    end
    methods
        function model = UncoupledReactionThermalModel()
            model.Reaction = ReactionModel();
            model.Thermal  = ThermalModel();
        end
    end
end
