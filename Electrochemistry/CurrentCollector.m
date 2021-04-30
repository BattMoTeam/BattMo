classdef CurrentCollector < ThermoElectronicComponent

    properties
        couplingTerm
    end
    
    methods
        
        function model = CurrentCollector(paramobj)
            
            model = model@ThermoElectronicComponent(paramobj);

            model = dispatchParams(model, paramobj, 'couplingTerm');
            
            % The parameter EffectiveElectronicConductivity in CurrentCollectorInputParams is given as scalar
            model.EffectiveElectronicConductivity = model.EffectiveElectronicConductivity*ones(model.G.cells.num, 1);
            
        end
        
        function state = updatejBcSource(model, state)
            
            state.jBcSource = state.jCoupling + state.jExternal;
            
        end
        
        function jExternal = setupExternalCoupling(model, phi, phiExternal);
            
            coupterm = model.couplingTerm;
            
            jExternal = phi*0.0; %NB hack to initialize zero ad
            
            sigmaeff = model.EffectiveElectronicConductivity;
            faces = coupterm.couplingfaces;
            bcval = phiExternal;
            [t, cells] = model.operators.harmFaceBC(sigmaeff, faces);
            jExternal(cells) = jExternal(cells) + t.*(bcval - phi(cells));
            
        end
        
    end
    
end
                             
