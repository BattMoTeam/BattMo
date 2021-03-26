classdef ActiveMaterial_ < CompositeModel
    
    methods
        
        function model = ActiveMaterial_(name)

            model = model@CompositeModel(name);
            names = {'T'       , ...
                     'SOC'     , ...
                     'phi'     , ... 
                     'cs'      , ...
                     'phiElyte', ...
                     'csElyte' , ...
                     'R'       , ...
                     'OCP'     , ...
                     'k'};
            model.names = names;
            model.vardims('cs') = 2;
            model.vardims('csElyte') = 2;
            
            % Alias: charge Carrier for cs{1}
            model = model.setAlias({'chargeCarrier', VarName({'.'}, 'cs', 2, 1)});
            model = model.setAlias({'chargeCarrierElyte', VarName({'.'}, 'csElyte', 2, 1)});
            
            fn = @ActiveMaterial.updateMaterialProperties;
            inputnames = {'chargeCarrier', 'T'};
            fnmodel = {'.'};
            
            propfunctions = model.propfunctions;
            
            propfunctions{end + 1} = PropFunction('OCP', fn, inputnames, fnmodel);
            propfunctions{end + 1} = PropFunction('k', fn, inputnames, fnmodel);
            
            fn = @ActiveMaterial.updateReactionRate;
            inputnames = {'T', 'phiElyte', 'phi', 'chargeCarrier', 'chargeCarrierElyte', 'OCP', 'k'};
            fnmodel = {'.'};
            propfunctions{end + 1} = PropFunction('R', fn, inputnames, fnmodel);

            model.propfunctions = propfunctions;

        end        
        
    end
    
end
