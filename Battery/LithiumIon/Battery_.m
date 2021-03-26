classdef Battery_ < CompositeModel

    methods

        function model = Battery_()

            model = model@CompositeModel('battery');
            
            names = {'T', ...
                     'SOC', ...
                     'prevstate', ...
                     'dt'};
            
            model.names = names;
            
            submodels = {};
            submodels{end + 1} = ElectroChemicalComponent_('elyte');
            submodels{end + 1} = Electrode_('ne');
            submodels{end + 1} = Electrode_('pe');
            
            model.SubModels = submodels;
            model.hasparent = false;

            model = model.initiateCompositeModel();
        end
        
    end
    
end
