classdef ActiveMaterial_ < CompositeModel
    
    methods
        
        function model = ActiveMaterial_(name)

            model = model@CompositeModel(name);
            names = {'T'             , ...
                     'SOC'           , ...
                     'phi'           , ... 
                     'cs'            , ...
                     'phiElectrolyte', ...
                     'csElectrolyte' , ...
                     'R'             , ...
                     'OCP'           , ...
                     'k'};
            model.names = names;
            model.vardims('cs') = 2;
            model.vardims('csElectrolyte') = 2;
            
            % Alias: charge Carrier for cs{1}
            model = model.setAlias({'chargeCarrier', VarName({'.'}, 'cs', 2, 1)});
            model = model.setAlias({'chargeCarrierElectrolyte', VarName({'.'}, 'csElectrolyte', 2, 1)});
            
            fn = @ActiveMaterial.updateMaterialProperties;
            inputnames = {'chargeCarrier', 'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('OCP', fn, inputnames, fnmodel);
            model = model.addPropFunction('k', fn, inputnames, fnmodel);
            
            fn = @ActiveMaterial.updateReactionRate;
            % inputnames = {'T', 'phiElectrolyte', 'phi', 'chargeCarrier', 'chargeCarrierElectrolyte', 'OCP', 'k'};
            inputnames = {'T','phi', 'phiElectrolyte', 'OCP', 'k'}; % for the moment, we do not consider dependance on concentration
            fnmodel = {'.'};
            model = model.addPropFunction('R', fn, inputnames, fnmodel);


        end        
        
    end
    
end
