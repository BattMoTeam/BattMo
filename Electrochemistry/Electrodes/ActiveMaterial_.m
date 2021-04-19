classdef ActiveMaterial_ < CompositeModel
    
    methods
        
        function model = ActiveMaterial_(name)

            model = model@CompositeModel(name);
            names = {'T'             , ...
                     'SOC'           , ...
                     'phiElectrode'           , ... 
                     'cElectrode'            , ...
                     'phiElectrolyte', ...
                     'cElectrolyte' , ...
                     'R'             , ...
                     'D'             , ...
                     'OCP'           , ...
                     'k'};
            model.names = names;
            
            % Alias: charge Carrier for cs{1}
            model = model.setAlias({'chargeCarrierElectrode', VarName({'.'}, 'cElectrode')});
            model = model.setAlias({'chargeCarrierElectrolyte', VarName({'.'}, 'cElectrolyte')});
            
            fn = @ActiveMaterial.updateDiffusionConductivityCoefficients;
            inputnames = {'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('D', fn, inputnames, fnmodel);
            model = model.addPropFunction('k', fn, inputnames, fnmodel);
            
            fn = @ActiveMaterial.updateOCP;
            inputnames = {'chargeCarrierElectrode', 'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('OCP', fn, inputnames, fnmodel);

            fn = @ActiveMaterial.updateReactionRate;
            % inputnames = {'T', 'phiElectrolyte', 'phiElectrode', 'chargeCarrierElectrode', 'chargeCarrierElectrolyte', 'OCP', 'k'};
            inputnames = {'T','phiElectrode', 'phiElectrolyte', 'OCP', 'k'}; % for the moment, we do not consider dependance on concentration
            fnmodel = {'.'};
            model = model.addPropFunction('R', fn, inputnames, fnmodel);

        end        
        
    end
    
end
