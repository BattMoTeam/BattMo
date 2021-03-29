classdef ElectronicComponent < PhysicalModel
    
    properties
        
        constants; % physical constants, see PhysicalConstants
        EffectiveElectronicConductivity;
    end

    methods
        
        function model = ElectronicComponent()
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

            flux   = state.j;
            bcflux = state.jBcSource;
            source = state.eSource;
            accum  = 0;
            
            chargeCons = assembleConservationEquation(model, flux, bcflux, source, accum);
            
            state.chargeCons = chargeCons;
            
        end
        
    end
end

