classdef Electrode_ < CompositeModel

    methods

        function model = Electrode_(name)
            
            model = model@CompositeModel(name);
            
            names = {model.names{:}, ...
                     'T' };
            
            model.SubModels{1} = ActiveElectroChemicalComponent_('main');
            model.SubModels{2} = ElectronicComponent_('cc');

            % Add potential coupling between the two
            main = model.getAssocModel('main');
            cc = model.getAssocModel('cc');
            
            inputnames = {VarName({'..', 'main'}, 'phi'), ...
                          VarName({'..', 'cc'}, 'phi')};
            fn = @Electrode.setupMainAndCCcoupling;
            fnmodel = {'..'};
            main = main.addPropFunction('jBcSource', fn, inputnames, fnmodel);
            cc = cc.addPropFunction('jBcSource', fn, inputnames, fnmodel);

            inputnames = {VarName({'..'}, 'T')};
            fn = @Electrode.updateT;
            fnmodel = {'..'};
            main = main.addPropFunction('T', fn, inputnames, fnmodel);
            cc = cc.addPropFunction('T', fn, inputnames, fnmodel);
            
            
            model = model.setSubModel('main', main);
            model = model.setSubModel('cc', cc);
            
        end
        
    end
    
end

