classdef NMC111Electrode_ < CompositeModel

    methods

        function model = NMC111Electrode_(name)

            model = model@CompositeModel(name);
            
            submodels = {};
            submodels{end + 1} = NMC111_('am');
            model.SubModels = submodels;
            
            model = model.setAlias({'cLi', VarName({'am'}, 'cLi')});
            
        end
        
    end
    
end

       