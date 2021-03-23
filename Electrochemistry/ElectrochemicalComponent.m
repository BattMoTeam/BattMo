classdef ElectrochemicalComponent < PhysicalModel
    
    properties
        
        constants 
        
        Volume
        EffectiveElectronicConductivity
        
    end

    
    methods
        
        function model = ElectrochemicalComponent()
            model = model@PhysicalModel([]);
            model.constants = PhysicalConstants();
        end

        function state = updateCurrent(model, state)
            
            sigmaeff = model.EffectiveElectronicConductivity;
            phi = state.phi;
            
            j = assembleFlux(model, phi, sigmaeff); 

            state.j = j;
            
        end
        
        function state = updateChargeConservation(model, state)
            
            state = model.updateCurrent(state);

            j = state.j;
            jBcSource = state.jBcSource;
            eSource   = state.eSource;
            
            massCons = assembleConservationEquation(model, j, jBcSource, eSource);
            
            state.massCons = massCons;
            
        end
        
    end
end

