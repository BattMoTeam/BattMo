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

        function [globalnames, localnames] = getModelPrimaryVarNames(model)
            localnames = {'phi'};
            globalnames = model.setupGlobalNames(localnames);
        end
        
        function [globalnames, localnames] = getModelVarNames(model)
            localnames = {'E', ...       % Potential at the end of collector
                          'j', ...       % Current density, [A/m2]
                          'OCP', ...     % Open-circuit potential [V];
                          'refOCP', ...  % Reference open circuit potential at standard temperature [V]
                         };
            globalnames = model.setupGlobalNames(localnames);
        end        

        function [globalnames, localnames] = getVarNames(model)
            [globalnames, localnames] = getVarNames@SimpleModel(model);
            globalnames{end + 1} = 'T';
            localnames{end + 1} = 'T';
        end

    end
end

