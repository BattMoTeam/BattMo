classdef ActiveElectroChemicalComponent < ElectroChemicalComponent
    
    properties
        
    end

    methods
        
        function model = ActiveElectroChemicalComponent()
            model = model@ElectroChemicalComponent([]);
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

