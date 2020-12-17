classdef currentCollector < SimpleModel
% Current collector model 
    
    properties
        E    % Potential at the end of collector
        eps
        sigma
        sigmaeff
    end
    
    methods
        function model = currentCollector(name, G, cells)
            
            model = model@SimpleModel(name);
            model.eps = 1;
            model.G = genSubGrid(G, cells);
            
        end
        
        function state = initializeState(model, state)
        % nothing to do here
        end

        function [namespaces, names] = getModelPrimaryVarNames(model)
            names = {'phi'};
            namespaces = model.assignCurrentNameSpace(names);
        end
        
        function [namespaces, names] = getModelVarNames(model)
            names = { 'j', ...      % Current density, [A/m2]
                     'OCP', ...     % Open-circuit potential [V];
                     'refOCP', ...  % Reference open circuit potential at standard temperature [V]
                    };
            namespaces = model.assignCurrentNameSpace(names);
        end        

        function [namespaces, names] = getVarNames(model)
            [namespaces, names] = getVarNames@SimpleModel(model);
            names = horzcat(names, {'T'});
            namespaces = horzcat(namespaces, {{}});
        end

    end
end

