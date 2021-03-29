classdef ActiveElectroChemicalComponent < ElectroChemicalComponent
    
    properties
        
    end

    methods
        
        function model = ActiveElectroChemicalComponent()
            model = model@ElectroChemicalComponent([]);
        end

        function state = updateIonAndCurrentSource(model, state)
            
            error('virtual function');
            
            % Uses reaction rate from active material to define the current and chargeCarrier sources:
            % Sketch of the method
            %
            % ionSourceName = model.ionSourceName;
            % R = state.ActiveMaterial.R;
            %
            % state.eSource = someFunctionOf(R);
            % state.(ionSourceName) = someFunctionOf(R);
            %
            
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

