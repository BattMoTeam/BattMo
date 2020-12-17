classdef currentCollector < SimpleModel
% Current collector model 
    
    properties
        eps
        sigma
        sigmaeff
    end
    
    methods
        function model = currentCollector(G, cells)
            
            model = model@SimpleModel();
            model.eps = 1;
            model.G = genSubGrid(G, cells);
            
        end

        function [namespaces, names] = getModelPrimaryVarNames(model)
            names = {'phi'};
            namespaces = model.assignCurrentNameSpace(names);
        end
        
        function [namespaces, names] = getModelVarNames(model)
            names = {'E', ...       % Potential at the end of collector
                     'j', ...       % Current density, [A/m2]
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

