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

            ne = model.getAssocModel('ne');
            pe = model.getAssocModel('pe');
            elyte = model.getAssocModel('elyte');
            
            fn = @Battery.UpdateTemp;
            inputnames = {VarName({'..'}, 'T')};
            fnmodel = {'..'};
            ne = ne.addPropFunction('T', fn, inputnames, fnmodel);
            pe = pe.addPropFunction('T', fn, inputnames, fnmodel);
            elyte = elyte.addPropFunction('T', fn, inputnames, fnmodel);
            
            model = model.setSubModel('ne', ne);
            model = model.setSubModel('pe', pe);
            model = model.setSubModel('elyte', elyte);            
            
            model = model.initiateCompositeModel();
        end
        
    end
    
end
