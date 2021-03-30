classdef ElectroChemicalComponent < ElectronicComponent
    
    properties

        % Names for book-keeping
        ionName
        ionFluxName 
        ionSourceName
        ionMassConsName
        ionAccumName
    end

    methods
        
        function model = ElectroChemicalComponent(paramobj)
            
            model = model@ElectronicComponent(paramobj);
            
            fdnames = {'ionName', ...
                       'ionFluxName', ...
                       'ionSourceName', ...
                       'ionMassConsName', ...
                       'ionAccumName'};
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
            
            ionName         = model.ionName;
            ionFluxName     = model.ionFluxName;
            ionSourceName   = model.ionSourceName;
            ionAccumName    = model.ionAccumName;
            ionMassConsName = model.ionMassConsName;
            
            flux   = state.(ionFluxName);
            source = state.(ionSourceName);
            accum  = state.(ionAccumName);
            bcflux = 0;
            
            masscons = assembleConservationEquation(model, flux, bcflux, source, accum);
            
            state.(ionMassConsName) = masscons;
            
        end
    end
end

