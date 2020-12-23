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

            model.pnames = {'phi'};
        
            % state variables
            names = {'phi', ...     % Potential
                     'j', ...       % Current density, [A/m2]
                     'OCP', ...     % Open-circuit potential [V];
                     'refOCP', ...  % Reference open circuit potential at standard temperature [V]
                    };
            model.names = names;
            model.aliases = {{'T', VarName({}, 'T')}};
            
        end
        
        function state = initializeState(model, state)
        % nothing to do here
        end
        
        function state = updateFlux(model, state)
            op = model.operators;
            sigmaeff = model.sigmaeff;
            
            phi = model.getProp(state, 'phi');
            
            j = - op.harmFace(sigmaeff).*op.Grad(phi); 
            
            state = model.setProp(state, 'j', j);
            
        end
        
    end
end
                             
