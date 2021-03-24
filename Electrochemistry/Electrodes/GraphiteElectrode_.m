classdef GraphiteElectrode_ < Electrode_

    methods

        function model = GraphiteElectrode_(name)

            model = model@Electrode_(name);

            submodels = {};
            submodels{end + 1} = Graphite_('am');
            model.SubModels = submodels;

            model = model.setAlias({'Li', VarName({'am'}, 'Li')});
           
        end
        
        
    end
    
end

       