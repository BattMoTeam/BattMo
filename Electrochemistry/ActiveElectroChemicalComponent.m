classdef ActiveElectroChemicalComponent < ElectroChemicalComponent
    
    properties
        
        ActiveMaterial
        
        volumeFraction
        porosity
        thickness
        
    end
    
    methods
        
        function model = ActiveElectroChemicalComponent(paramobj)
        % shortcut used here:
        % am = ActiveMaterial
            
            model = model@ElectroChemicalComponent(paramobj);
            
            % Setup ActiveMaterial component
            paramobj.am.G = model.G;
            
            model.ActiveMaterial = setupActiveMaterial(paramobj.am);

            % setup volumeFraction, porosity, thickness
            nc = model.G.cells.num;
            volumeFraction = model.ActiveMaterial.volumeFraction*ones(nc, 1);
            model.volumeFraction = volumeFraction;
            model.porosity = 1 - model.volumeFraction;
            model.thickness = 10e-6;
            
            % setup effective electronic conductivity
            econd = model.ActiveMaterial.electronicConductivity;
            model.EffectiveElectronicConductivity = econd .* volumeFraction.^1.5;
            
        end

        function state = updateIonAndCurrentSource(model, state)
            
            ccSourceName = model.chargeCarrierSourceName;
            
            R = state.ActiveMaterial.R;
            
            state.eSource = R;
            state.(ccSourceName) = -R;
            
        end
        
        function state = updateDiffusionCoefficient(model, state)
            
            D = state.ActiveMaterial.D;
            state.D = D .* model.volumeFraction .^1.5;
            
        end
        
        function state = updateChargeCarrier(model, state)
            
            state.c = state.ActiveMaterial.c;
            
        end 
        
        function state = updatePhi(model, state)
            state.phi = state.ActiveMaterial.phi;
        end         
        
        function state = updateT(model, state)
            state.ActiveMaterial.T = state.T;
        end

        function state = updatejBcSource(model, state)
            state.jBcSource = state.jCoupling;
        end
        
    end
end

