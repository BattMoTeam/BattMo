classdef Electrode_ < CompositeModel

    methods

        function model = Electrode_(name)
            
            model = model@CompositeModel(name);
            
            names = {model.names{:}, ...
                     'T' };
            
            model.SubModels{1} = ActiveElectroChemicalComponent_('aecm');
            model.SubModels{2} = ElectronicComponent_('cc');

            % Add potential coupling between the two
            aecm = model.getAssocModel('aecm');
            cc = model.getAssocModel('cc');
            
            inputnames = {VarName({'..', 'aecm'}, 'phi'), ...
                          VarName({'..', 'cc'}, 'phi')};
            fn = @Electrode.setupCCcoupling;
            fnmodel = {'..'};
            aecm = aecm.addPropFunction('jBcSource', fn, inputnames, fnmodel);
            cc = cc.addPropFunction('jBcSource', fn, inputnames, fnmodel);

            inputnames = {VarName({'..'}, 'T')};
            fn = @Electrode.updateT;
            fnmodel = {'..'};
            aecm = aecm.addPropFunction('T', fn, inputnames, fnmodel);
            cc = cc.addPropFunction('T', fn, inputnames, fnmodel);
            
            
            model = model.setSubModel('aecm', aecm);
            model = model.setSubModel('cc', cc);
            
        end
        
    end
    
end

