classdef Electrode < PhysicalModel
    
    properties
        
        % submodels
        ElectrodeActiveComponent % is a ActiveElectroChemicalComponent
        CurrentCollector % is a ElectronicComponent
        
    end

    methods
        
        function model = Electrode()
            model = PhysicalModel([]);
            
            % Setup the two components
            
        end
        
        function state  = setupCurrentCollectorCoupling(model, state)
        % setup electrical coupling between current collector and the electrode active component
            
            error('to be implemented')
            
            eac_phi = state.ElectrodeActiveComponent.phi;
            cc_phi = state.CurrentCollector.phi;
            
            state.ElectrodeActiveComponent.jBcSource = someFunctionOf(eac_phi, cc_phi);
            state.CurrentCollector.jBcSource = someFunctionOf(eac_phi, cc_phi);
            
        end
        
        function state = updateT(model, state)
            state.ElectrodeActiveComponent.T = state.T;
            state.CurrentCollector.T = state.T;
        end        
        
    end    
end

