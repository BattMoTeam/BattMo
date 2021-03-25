classdef ActiveMaterial_ < ComponentModel
    
    methods
        
        function model = ActiveMaterial_(name)

            model = model@ComponentModel(name);
            names = {'T'    , ...
                     'SOC'  , ...
                     'cs'   , ...
                     'theta', ...
                     'D'    , ...
                     'OCP'  , ...
                     'k'    , ...
                     'phi'};
            model.names = names;

            % Alias: Li for cs{1}
            model = model.setAlias({'Li', VarName({'.'}, 'cs', 2, 1)});
            
            fn = @ActiveMaterial.updateMaterialProperties;
            inputnames = {'c', 'T'};
            fnmodel = {'.'};
            
            propfunctions = model.propfunctions;
            
            propfunctions{end + 1} = PropFunction('theta', fn, inputnames, fnmodel);
            propfunctions{end + 1} = PropFunction('D', fn, inputnames, fnmodel);
            propfunctions{end + 1} = PropFunction('OCP', fn, inputnames, fnmodel);
            propfunctions{end + 1} = PropFunction('k', fn, inputnames, fnmodel);
            
            model.propfunctions = propfunctions;

        end        
        
    end
    
end
