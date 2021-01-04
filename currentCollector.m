classdef currentCollector < ComponentModel
% Current collector model 
    
    properties
        con = physicalConstants();
        
        E    % Potential at the end of collector
        eps
        sigma
        sigmaeff
    end
    
    methods
        function model = currentCollector(name, G, cells)
            
            model = model@ComponentModel(name);
            model.G = genSubGrid(G, cells);

            model.sigma = 100;
            model.eps = ones(model.G.cells.num, 1);
            model.sigmaeff = (model.sigma) .* (model.eps).^1.5;            

            % setup operators
            model.operators = localSetupOperators(model.G);
            
            model.pnames = {'phi'};
        
            % state variables
            names = {'phi', ...     % Potential
                     'T', ...       % Temperature
                     'j', ...       % Current density, [A/m2]
                     'jBcSource', ...
                     'chargeCont', ...
                     'OCP', ...     % Open-circuit potential [V];
                     'refOCP', ...  % Reference open circuit potential at standard temperature [V]
                    };
            model.names = names;
            
            propfunctions = {};
            
            % setup updating function for j
            name = 'j';
            updatefn = @(model, state) model.updateFlux(state);
            propfunction = PropFunction(name, updatefn, '.');
            propfunctions{end + 1} = propfunction;
            
            % setup update property function for charge continuity (chargeCont)
            name = 'chargeCont';
            updatefn = @(model, state) model.updateChargeCont(state);
            propfunction = PropFunction(name, updatefn, '.');
            propfunctions{end + 1} = propfunction;
            
            
            model.propfunctions = propfunctions;
            
        end
        
        function state = initializeState(model, state)
            % Used only in debugging for the moment
            state = model.initiateState(state);
        end
        
        function state = updateFlux(model, state)
            op = model.operators;
            sigmaeff = model.sigmaeff;
            
            [phi, state] = model.getUpdatedProp(state, 'phi');
            
            j = - op.harmFace(sigmaeff).*op.Grad(phi); 
            
            state = model.setProp(state, 'j', j);
            
        end
        
        function state = updateChargeCont(model, state)
            
            op = model.operators;
            
            [j, state]         = model.getUpdatedProp(state, 'j');
            [jBcSource, state] = model.getUpdatedProp(state, 'jBcSource');
            
            chargeCont = (op.Div(j) - jBcSource)./ model.G.cells.volumes./model.con.F;
            
            state = model.setProp(state, 'chargeCont', chargeCont);
            
        end
         
    end
end
                             
