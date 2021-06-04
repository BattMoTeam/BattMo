classdef ElectroChemicalComponent < ElectronicComponent
    
    properties

        % Names for book-keeping
        chargeCarrierName
        chargeCarrierFluxName 
        chargeCarrierSourceName
        chargeCarrierMassConsName
        chargeCarrierAccumName
        
        EffectiveDiffusionCoefficient
        
    end

    methods
        
        function model = ElectroChemicalComponent(paramobj)
            
            model = model@ElectronicComponent(paramobj);
            
            fdnames = {'chargeCarrierName', ...
                       'chargeCarrierFluxName', ...
                       'chargeCarrierSourceName', ...
                       'chargeCarrierMassConsName', ...
                       'chargeCarrierAccumName'};
            model = dispatchParams(model, paramobj, fdnames);
            
        end

        function state = updateChargeCarrierFlux(model, state)
            
            ccFluxName = model.chargeCarrierFluxName;

            D = model.EffectiveDiffusionCoefficient;
            
            c = state.c;

            ccflux = assembleFlux(model, c, D);
            
            state.(ccFluxName) = ccflux;
            
        end
        
        function state = updateMassConservation(model, state)
            
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

