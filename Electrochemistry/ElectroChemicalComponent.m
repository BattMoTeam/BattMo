classdef ElectroChemicalComponent < ElectronicComponent
    
    properties

        % Names for book-keeping
        chargeCarrierName
        chargeCarrierFluxName 
        chargeCarrierSourceName
        chargeCarrierMassConsName
        chargeCarrierAccumName
        
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

        function state = updateDiffusionCoefficient(model, state)
            error('virtual function : diffusion coefficient is model dependent')
        end

        function state = updateChargeCarrierFlux(model, state)
            
            ccFluxName = model.chargeCarrierFluxName;
            
            D = state.D;
            c = state.c;

            ccflux = assembleFlux(model, c, D);
            
            %% apply scaling (maybe not the right place but consistent with assembleConservationEquation - at
            %% least when this comment was written...)
            F = model.constants.F;
            ccflux = ccflux*F;
            
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
            bcflux = 0;
            
            masscons = assembleConservationEquation(model, flux, bcflux, source, accum);
            
            state.(ccMassConsName) = masscons;
            
        end
    end
end

