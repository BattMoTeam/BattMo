classdef ActiveElectroChemicalComponent < ElectroChemicalComponent
    
    properties
        ActiveMaterial
    end
    
    methods
        
        function model = ActiveElectroChemicalComponent(paramobj)
            model = model@ElectroChemicalComponent(paramobj);
            paramobj.am.G = model.G;
            model.ActiveMaterial = ActiveMaterial(paramobj.am);
        end

        function state = updateIonAndCurrentSource(model, state)
            
            R = state.ActiveMaterial.R;
            
            state.eSource = R;
            state.(ionSourceName) = -R;
            
        end
        
        function state = updateChargeCarrier(model, state)
            state.ActiveMaterial.(ionName) = state.(ionName);
        end 
        
        function state = updatePhi(model, state)
            state.ActiveMaterial.phi = state.phi;            
        end         
        
        function state = updateT(model, state)
            state.ActiveMaterial.T = state.T;
        end
        
        
    end
end

