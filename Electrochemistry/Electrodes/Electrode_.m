classdef Electrode_ < CompositeModel

    methods

        function model = Electrode_(name)
            
            model = model@CompositeModel(name);
            
            names = {model.names{:}, ...
                     'T' };
            
            model.SubModels{1} = ElectrodeActiveComponent_('eac');
            model.SubModels{2} = CurrentCollector_('cc');

            % Add coupling between the two components
            inputnames = {VarName({'..', 'eac'}, 'phi'), ...
                          VarName({'..', 'cc'}, 'phi')};
            fn = @Electrode.setupCoupling;
            fnmodel = {'..'};
            model = model.addPropFunction({'eac', 'jCoupling'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'cc', 'jCoupling'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'cc', 'eSource'}, fn, inputnames, fnmodel);
            
            % dispatch temperature
            inputnames = {VarName({'..'}, 'T')};
            fn = @Electrode.updateT;
            fnmodel = {'..'};
            model = model.addPropFunction({'eac', 'T'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'cc', 'T'}, fn, inputnames, fnmodel);
            
        end
        
    end
    
end

