classdef ReactionModel2 < ReactionModel

    methods

        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@ReactionModel(model);
            
            varnames = {};
            % potential in electrode
            varnames{end + 1} = 'T';
            
            model = model.registerVarNames(varnames);
            
            fn = @ReactionModel.updateOCP;
            inputnames = {'c_s', 'T'};
            model = model.registerPropFunction({'OCP', fn, inputnames});
            
        end

    end
    
end
