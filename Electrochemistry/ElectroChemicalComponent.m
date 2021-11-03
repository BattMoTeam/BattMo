classdef ElectroChemicalComponent < ElectronicComponent
    
    properties

        EffectiveDiffusionCoefficient % Effective diffusion coefficient
        
        
    end

    methods
        
        function model = ElectroChemicalComponent(paramobj)
        % Here, :code:`paramobj` is instance of :class:`ElectronicChemicalInputParams <Electrochemistry.ElectronicChemicalInputParams>`
        
            model = model@ElectronicComponent(paramobj);
            
        end

        function state = updateMassFlux(model, state)
        % Assemble diffusion flux which is stored in :code:`state.Flux`

            D = model.EffectiveDiffusionCoefficient;
            
            c = state.c;

            massflux = assembleFlux(model, c, D);
            
            state.massFlux = massflux;
            
        end
        
        function state = updateMassConservation(model, state)
        % Assemble residual of the mass conservation equation which is stored in :code:`state.massCons`
            
            flux   = state.massFlux;
            source = state.massSource;
            accum  = state.massAccum;
            bcsource = 0;
            
            masscons = assembleConservationEquation(model, flux, bcsource, source, accum);
            
            state.massCons = masscons;
            
        end
        
    end
end

