classdef Battery_ < CompositeModel

    methods

        function model = Battery_()

            model = model@CompositeModel('battery');
            
            names = {'SOC', ...
                     'prevstate', ...
                     'energyCons', ...
                     'dt'};
            
            model.names = names;
            
            submodels = {};
            submodels{end + 1} = Electrolyte_('elyte');
            submodels{end + 1} = Electrode_('ne');
            submodels{end + 1} = Electrode_('pe');
            submodels{end + 1} = ThermalComponent_('thermal');
            
            model.SubModels = submodels;
            
            %% update temperatures (dispatching)
            fn = @Battery.updateTemperature;
            
            fnmodel = {'..', '..'};
            inputnames = {VarName({'..', '..', 'thermal'}, 'T')};
            model = model.addPropFunction({'ne', 'eac', 'T'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'ne', 'cc', 'T'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'pe', 'eac', 'T'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'pe', 'cc', 'T'}, fn, inputnames, fnmodel);            
            fnmodel = {'..'};
            inputnames = {VarName({'..', 'thermal'}, 'T')};            
            model = model.addPropFunction({'elyte', 'T'}, fn, inputnames, fnmodel);
           
                  
            %% setup couplings
            
            fn = @Battery.updateElectrodeCoupling;
            
            clear inputnames;
            inputnames{1} = VarName({'elyte'}, 'chargeCarrier');
            inputnames{1}.isNamingRelative = false;
            inputnames{2} = VarName({'elyte'}, 'phi');
            inputnames{2}.isNamingRelative = false;
            
            fnmodel = {'..', '..', '..'};
            model = model.addPropFunction({'ne', 'eac', 'am', 'phiElectrolyte'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'ne', 'eac', 'am', 'chargeCarrierElectrolyte'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'pe', 'eac', 'am', 'phiElectrolyte'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'pe', 'eac', 'am', 'chargeCarrierElectrolyte'}, fn, inputnames, fnmodel);
            
            fn = @Battery.updateElectrolyteCoupling;
            
            clear inputnames;
            inputnames{1} = VarName({'ne', 'eac', 'am'}, 'R');
            inputnames{1}.isNamingRelative = false;
            inputnames{2} = VarName({'pe', 'eac', 'am'}, 'R');
            inputnames{2}.isNamingRelative = false;
            
            fnmodel = {'..'};
            model = model.addPropFunction({'elyte', 'chargeCarrierSource'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'elyte', 'eSource'}, fn, inputnames, fnmodel);
            
                
            fn = @Battery.updateAccumTerms;
            
            inputnames = {VarName({'.'}, 'chargeCarrier'), ...
                          VarName({'..', 'ne', 'eac', 'am'}, 'chargeCarrier'), ...
                          VarName({'..', 'pe', 'eac', 'am'}, 'chargeCarrier')};
            
            fnmodel = {'..'};
            model = model.addPropFunction({'elyte', 'chargeCarrierAccum'}, fn, inputnames, fnmodel);            


            inputnames = {VarName({'..', '..', 'elyte'}, 'chargeCarrier'), ...
                          VarName({'am'}, 'chargeCarrier'), ...
                          VarName({'..', '..', 'pe', 'eac', 'am'}, 'chargeCarrier')};
            
            fnmodel = {'..', '..'};
            
            model = model.addPropFunction({'ne', 'eac', 'chargeCarrierAccum'}, fn, inputnames, fnmodel);
            
            inputnames = {VarName({'..', '..', 'elyte'}, 'chargeCarrier'), ...
                          VarName({'..', '..', 'ne', 'eac', 'am'}, 'chargeCarrier'), ...
                          VarName({'am'}, 'chargeCarrier')};
            
            fnmodel = {'..', '..'};
            
            model = model.addPropFunction({'pe', 'eac', 'chargeCarrierAccum'}, fn, inputnames, fnmodel);                        
            
            %% update Thermal accumulation terms
            
            fn = @Battery.updateThermalAccumTerms;
            fnmodel = {'..'};
            inputnames = {VarName({'thermal'}, 'T')};
            model = model.addPropFunction({'thermal', 'accumHeat'}, fn, inputnames, fnmodel);
            
            %% update Thermal Ohmic Terms
            
            fn = @Battery.updateThermalOhmicSourceTerms;
            fnmodel = {'..'};
            inputnames = {VarName({'elyte'}, 'j')     , ...
                          VarName({'ne', 'cc'}, 'j')  , ...
                          VarName({'ne', 'eac'}, 'j') , ...
                          VarName({'pe', 'cc'}, 'j')  , ...
                          VarName({'pe', 'eac'}, 'j')};
            model = model.addPropFunction({'thermal', 'jHeatOhmSource'}, fn, inputnames, fnmodel);
            model = model.addPropFunction({'thermal', 'jHeatBcSource'}, fn, inputnames, fnmodel);
            
            fn = @Battery.updateThermalSourceTerms;
            fnmodel = {'..'};
            inputnames = {VarName({'thermal'}, 'jHeatSource')};
            model = model.addPropFunction({'thermal', 'jHeatSource'}, fn, inputnames, fnmodel);
                          
                          
            %% setup external coupling at positive and negative electrodes
            
            fn = @Battery.setupExternalCouplingNegativeElectrode;
            inputnames = {'phi'};
            fnmodel = {'.'};
                      
            model = model.addPropFunction({'ne', 'cc', 'jExternal'}, fn, inputnames, fnmodel);
            
            fn = @Battery.setupExternalCouplingPositiveElectrode;
            inputnames = {'phi', ...
                          'E'};
            fnmodel = {'.'};
            model = model.addPropFunction({'pe', 'cc', 'jExternal'}, fn, inputnames, fnmodel);
            
            model = model.initiateCompositeModel();
            
        end
        
    end
    
end
