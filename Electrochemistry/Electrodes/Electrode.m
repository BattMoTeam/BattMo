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
        
    end

    methods
       
        function state = updateCurrent(model, state)
            
            sigmaeff = model.EffectiveElectronicConductivity;
            phi = state.ActiveMaterial.phi;
            
            j = assembleFlux(model, phi, sigmaeff); 

            state.j = j;
            
        end

        function state = updateIonFlux(model, state)
        % update Ion fluxes in electrode (this depends on the type of electrode - because it depends on the name of the ion!)
            error('virtual function')
        end
        
    end    
end

