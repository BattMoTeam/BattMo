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
            submodels{end + 1} = Electrolyte_('elyte');
            submodels{end + 1} = Electrode_('ne');
            submodels{end + 1} = Electrode_('pe');
            
            model.SubModels = submodels;
            
            fn = @Battery.UpdateT;
            inputnames = {VarName({'..'}, 'T')};
            fnmodel = {'..'};
            model = model.addPropFunction({'ne', 'T'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'pe', 'T'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'elyte', 'T'}, fn, inputnames, fnmodel);
            
            fn = @Battery.setupElectrodeCoupling;
            
            clear inputnames;
            inputnames{1} = VarName({'elyte'}, 'chargeCarrier');
            inputnames{1}.isNamingRelative = false;
            inputnames{2} = VarName({'elyte'}, 'phi');
            inputnames{2}.isNamingRelative = false;
            
            fnmodel = {'..', '..', '..'};
            model = model.addPropFunction({'ne', 'eac', 'am', 'phiElyte'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'ne', 'eac', 'am', 'chargeCarrierElyte'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'pe', 'eac', 'am', 'phiElyte'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'pe', 'eac', 'am', 'chargeCarrierElyte'}, fn, inputnames, fnmodel);
            
            fn = @Battery.setupElectrolyteCoupling;
            
            clear inputnames;
            inputnames{1} = VarName({'ne', 'eac', 'am'}, 'R');
            inputnames{1}.isNamingRelative = false;
            inputnames{1} = VarName({'pe', 'eac', 'am'}, 'R');
            inputnames{2}.isNamingRelative = false;
            
            fnmodel = {'..'};
            model = model.addPropFunction({'elyte', 'chargeCarrierSource'}, fn, inputnames, fnmodel);
            
            model = model.initiateCompositeModel();
            
            
        end
        
    end
    
end
