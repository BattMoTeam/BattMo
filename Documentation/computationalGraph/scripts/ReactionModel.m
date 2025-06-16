classdef ReactionModel < BaseModel

    methods

        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            % potential in electrode
            varnames{end + 1} = 'phi_s';
            % charge carrier concentration in electrode - value at surface
            varnames{end + 1} = 'c_s';
            % potential in electrolyte
            varnames{end + 1} = 'phi_e';
            % charge carrier concentration in electrolyte
            varnames{end + 1} = 'c_e';
            % Electrode over potential
            varnames{end + 1} = 'eta';
            % Intercalation flux [mol s^-1 m^-2]
            varnames{end + 1} = 'R';
            % OCP [V]
            varnames{end + 1} = 'OCP';
            % Reaction rate coefficient [A m^-2]
            varnames{end + 1} = 'j';
            
            model = model.registerVarNames(varnames);
            
            fn = @ReactionModel.updateReactionRateCoefficient;
            inputnames = {'c_e', 'c_s'};
            model = model.registerPropFunction({'j', fn, inputnames});

            fn = @ReactionModel.updateOCP;
            inputnames = {'c_s'};
            model = model.registerPropFunction({'OCP', fn, inputnames});
            
            fn = @ReactionModel.updateEta;
            inputnames = {'phi_e', 'phi_s', 'OCP'};            
            model = model.registerPropFunction({'eta', fn, inputnames});

            fn = @ReactionModel.updateReactionRate;
            inputnames = {'eta', 'j'};
            model = model.registerPropFunction({'R', fn, inputnames});
            
        end

    end
    
end
