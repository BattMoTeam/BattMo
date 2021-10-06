classdef ElectroChemicalComponent < ElectronicComponent
    
    properties

        EffectiveDiffusionCoefficient % Effective diffusion coefficient

        chargeCarrierName
        chargeCarrierFluxName 
        chargeCarrierSourceName
        chargeCarrierMassConsName
        chargeCarrierAccumName
        
        
    end

    methods
        
        function model = ElectroChemicalComponent(paramobj)
        % Here, :code:`paramobj` is instance of :class:`ElectronicChemicalInputParams <Electrochemistry.ElectronicChemicalInputParams>`
        
            model = model@ElectronicComponent(paramobj);
            
            fdnames = {'chargeCarrierName'};
            model = dispatchParams(model, paramobj, fdnames);
            
            ccname = model.chargeCarrierName;
            model.chargeCarrierFluxName     = sprintf('%sFlux', ccname);
            model.chargeCarrierSourceName   = sprintf('%sSource', ccname);
            model.chargeCarrierMassConsName = 'massCons';
            model.chargeCarrierAccumName    = sprintf('%sAccum', ccname);
            
        end

        function state = updateChargeCarrierFlux(model, state)
        % Assemble diffusion flux which is stored in :code:`state.Flux`
            ccFluxName = model.chargeCarrierFluxName;

            D = model.EffectiveDiffusionCoefficient;
            
            c = state.c;

            ccflux = assembleFlux(model, c, D);
            
            state.(ccFluxName) = ccflux;
            
        end
        
        function state = updateMassConservation(model, state)
        % Assemble residual of the mass conservation equation which is stored in :code:`state.massCons`
            
            ccName         = model.chargeCarrierName;
            ccFluxName     = model.chargeCarrierFluxName;
            ccSourceName   = model.chargeCarrierSourceName;
            ccAccumName    = model.chargeCarrierAccumName;
            ccMassConsName = model.chargeCarrierMassConsName;
            
            flux   = state.(ccFluxName);
            source = state.(ccSourceName);
            accum  = state.(ccAccumName);
            bcsource = 0;
            
            masscons = assembleConservationEquation(model, flux, bcsource, source, accum);
            
            state.(ccMassConsName) = masscons;
            
        end
    end
end

