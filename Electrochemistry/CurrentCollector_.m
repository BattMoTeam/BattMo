classdef CurrentCollector_ < ElectronicComponent_
    
    methods

        function model = CurrentCollector_(name)

            model = model@ElectronicComponent_(name);
            
            names = {model.names{:}, ...
                     'jCoupling', ...
                     'jExternal'};
            model.names = names;
            
            fn = @CurrentCollector.updatejBcSource;
            inputnames = {'jCoupling', 'jExternal'};
            fnmodel = {'.'};
            model = model.addPropFunction('jBcSource', fn, inputnames, fnmodel);            
            
        end
        
    end
    
end

       