classdef NMC111Electrode_ < Electrode_

    methods

        function model = NMC111Electrode_(name)

            model = model@Electrode_(name);
            
            submodels = {};
            submodels{end + 1} = NMC111_('am');
            model.SubModels = submodels;
            
            model = model.setAlias({'Li', VarName({'am'}, 'Li')});
            
        end
        
    end
    
end

       