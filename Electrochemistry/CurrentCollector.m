classdef CurrentCollector < ElectronicComponent

    properties
        
        couplingTerm

        thermalConductivity
        heatCapacity
        
    end
    
    methods
        
        function model = CurrentCollector(paramobj)
            
            model = model@ElectronicComponent(paramobj);

            fdnames = {'couplingTerm', ...
                       'thermalConductivity', ...
                       'heatCapacity'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
            % The parameter EffectiveElectricalConductivity in CurrentCollectorInputParams is given as scalar
            model.EffectiveElectricalConductivity = model.EffectiveElectricalConductivity*ones(model.G.cells.num, 1);
            
        end
        
        function state = updatejBcSource(model, state)
            
            state.jBcSource = state.jCoupling + state.jExternal;
            
        end
        
        function jExternal = setupExternalCoupling(model, phi, phiExternal);
            
            coupterm = model.couplingTerm;
            
            jExternal = phi*0.0; %NB hack to initialize zero ad
            
            sigmaeff = model.EffectiveElectricalConductivity;
            faces = coupterm.couplingfaces;
            bcval = phiExternal;
            [t, cells] = model.operators.harmFaceBC(sigmaeff, faces);
            jExternal(cells) = jExternal(cells) + t.*(bcval - phi(cells));
            
        end
        
    end
    
end
                             
