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
            
            ionName = model.ionName;
            ionFluxName = model.ionFluxName;
            
            D = state.D;
            c = state.(ionName);

            ionflux = assembleFlux(model, c, D);
            
            F = model.constants.F;
            ionflux = ionflux*F;
            
            state.(ionFluxName) = ionflux;
        end

        function state = updateCurrent(model, state)
        % typically will be overloaded to take into account charge Carrier presence
            state = updateCurrent@ElectronicComponent(model, state)
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

