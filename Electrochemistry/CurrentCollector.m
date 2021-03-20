classdef CurrentCollector < ComponentModel
% Current collector model 
    
    properties
        constants = PhysicalConstants();
        
        E    % Potential at the end of collector
        volumeFraction
        electronicConductivity
        effectiveElectronicConductivity
    end
    
    methods
        function model = CurrentCollector(name, G, cells)
            
            model = model@ComponentModel(name);
            model.G = genSubGrid(G, cells);

            model.electronicConductivity = 100;
            model.volumeFraction = ones(model.G.cells.num, 1);
            model.effectiveElectronicConductivity = (model.electronicConductivity) .* (model.volumeFraction).^1.5;            

            % setup operators
            model.operators = localSetupOperators(model.G);
            
            model.pnames = {'phi'};
        
            % state variables
            names = {'phi', ...     % Potential
                     'T', ...       % Temperature
                     'jBcSource', ...
                     'chargeCont', ...
                    };
            model.names = names;
            model = model.setupVarDims();
            
            % setup update property function for charge continuity (chargeCont)
            name = 'chargeCont';
            updatefn = @(model, state) model.updateChargeCont(state);
            inputnames = {'phi', 'jBcSource'};
            model = model.addPropFunction(name, updatefn, inputnames, {'.'});
            
        end
        
        function state = updateChargeCont(model, state)
            
            op = model.operators;
            sigmaeff = model.effectiveElectronicConductivity;
         
            phi = state.phi;
            jBcSource = state.jBcSource;
            
            j = - op.harmFace(sigmaeff).*op.Grad(phi); 
            
            chargeCont = (op.Div(j) - jBcSource)./ model.G.cells.volumes./model.constants.F;
            
            state.chargeCont = chargeCont;
            
        end
         
    end
end
                             
