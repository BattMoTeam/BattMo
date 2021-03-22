classdef ElectrochemicalComponent < PhysicalModel
    
    properties
        Volume
        EffectiveElectronicConductivity
    end

    function state = updateChargeConservation(model, state)
        
        phi       = state.phi;
        jBcSource = state.jBcSource;
        eSource   = state.eSource;
        
        massCons = assembleChargeConservationEquation(model, phi, jBcSource, eSource);
        
        state.massCons = massCons;
        
    end
end

