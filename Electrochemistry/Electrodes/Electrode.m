classdef Electrode < ElectrochemicalComponent
    
    properties
        
        % Design properties
        thickness       % Thickness,        [m]
        volumeFraction  % Volume fraction,  [-]
        porosity        % Porosity,         [-]
        mass            % Mass,             [kg]
        volume          % Volume,           [m3]
        area            % Area,             [m2]
        
        % SubModels
        ActiveMaterial
        Binder              % Binder object
        ConductingAdditive  % Conducting additive object        
        Electrolyte         % Electrolyte data structure

        % names for book-keeping
        ionName
        ionFluxName 
        ionSourceName
        ionMassConsName
        ionAccumName
    end

    methods
       
        function state = updateCurrent(model, state)
            
            sigmaeff = model.EffectiveElectronicConductivity;
            phi = state.ActiveMaterial.phi;
            
            j = assembleFlux(model, phi, sigmaeff); 

            state.j = j;
            
        end

        function state = updateIonFlux(model, state)
            ionName = model.ionName;
            ionFluxName = model.ionFluxName;
            
            D = state.ActiveMaterial.D;
            c = state.ActiveMaterial.(ionName);
            Deff = D .* model.volumeFraction .^1.5;
            
            ionflux = assembleFlux(model, c, Deff);
            
            state.(ionFluxName) = ionflux;
        end
  
        function state = updateMassConservationEquation(model, state)
            
            ionName       = model.ionName;
            ionFluxName   = model.ionFluxName;
            ionSourceName = model.ionSourceName;
            
            flux   = state.(ionFluxName);
            source = state.(ionSourceName);
            accum  = state.(ionAccumName);
            bcflux = 0;
            
            masscons = assembleConservationEquation(model, flux, bcflux, source, accum);
            
            state.(ionMassConsName) = masscons;
            
        end
        
    end    
end

