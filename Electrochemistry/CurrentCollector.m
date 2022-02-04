classdef CurrentCollector < ElectronicComponent

    properties
        
        couplingTerm

        thermalConductivity
        heatCapacity
        density
    end
    
    methods
        
        function model = CurrentCollector(paramobj)
            
            model = model@ElectronicComponent(paramobj);

            fdnames = {'couplingTerm'        , ...
                       'thermalConductivity' , ...
                       'heatCapacity'        , ...
                       'density'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
            % The parameter EffectiveElectricalConductivity in CurrentCollectorInputParams is given as scalar
            model.EffectiveElectricalConductivity = model.EffectiveElectricalConductivity*ones(model.G.cells.num, 1);
            
        end
        
        function state = updatejBcSource(model, state)
            
            state.jBcSource = state.jCoupling + state.jExternal;
            state.jFaceBc = state.jFaceCoupling + state.jFaceExternal;
            
        end
        
        function [jExternal, jFaceExternal] = setupExternalCoupling(model, phi, phiExternal);
            
            coupterm = model.couplingTerm;
            
            jExternal = phi*0.0; %NB hack to initialize zero ad
            
            sigmaeff = model.EffectiveElectricalConductivity;
            faces = coupterm.couplingfaces;
            bcval = phiExternal;
            [t, cells] = model.operators.harmFaceBC(sigmaeff, faces);
            current = t.*(bcval - phi(cells));
            jExternal(cells) = jExternal(cells) + current;
            
            G = model.G;
            nf = G.faces.num;
            %sgn = model.operators.sgn;
            zeroFaceAD = model.AutoDiffBackend.convertToAD(zeros(nf, 1), phi);
            jFaceExternal = zeroFaceAD;
            sgn = 2*(cells == G.faces.neighbors(faces,1))-1;
            jFaceExternal(faces) = -sgn.*current;
            
        end
        
    end
    
end
                             
