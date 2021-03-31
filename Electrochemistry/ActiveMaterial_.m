classdef ActiveMaterial_ < CompositeModel
    
    methods
        
        function model = ActiveMaterial_(name)

            model = model@CompositeModel(name);
            names = {'T'             , ...
                     'SOC'           , ...
                     'phi'           , ... 
                     'c'            , ...
                     'phiElectrolyte', ...
                     'cElectrolyte' , ...
                     'R'             , ...
                     'D'             , ...
                     'OCP'           , ...
                     'k'};
            model.names = names;
            
            % Alias: charge Carrier for cs{1}
            model = model.setAlias({'chargeCarrier', VarName({'.'}, 'c')});
            model = model.setAlias({'chargeCarrierElectrolyte', VarName({'.'}, 'cElectrolyte')});
            
            fn = @ActiveMaterial.updateMaterialProperties;
            inputnames = {'chargeCarrier', 'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('OCP', fn, inputnames, fnmodel);
            model = model.addPropFunction('D', fn, inputnames, fnmodel);
            model = model.addPropFunction('k', fn, inputnames, fnmodel);
            
            fn = @ActiveMaterial.updateReactionRate;
            % inputnames = {'T', 'phiElectrolyte', 'phi', 'chargeCarrier', 'chargeCarrierElectrolyte', 'OCP', 'k'};
            inputnames = {'T','phi', 'phiElectrolyte', 'OCP', 'k'}; % for the moment, we do not consider dependance on concentration
            fnmodel = {'.'};
            model = model.addPropFunction('R', fn, inputnames, fnmodel);


        end        
        
    end
    
end
