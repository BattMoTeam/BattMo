classdef ThermoElectroChemicalComponent < ElectroChemicalComponent
    
    properties
        
        EffectiveThermalConductivity
        heatCapacity % in [J][K]^-1[m]^-3
    
    end

    methods
        
        function model = ThermoElectroChemicalComponent(paramobj)
            
            model = model@ElectroChemicalComponent(paramobj);
            
            fdnames = {'EffectiveThermalConductivity', ...
                       'heatCapacity'};
            model = dispatchParams(model, paramobj, fdnames);
            
            nc = model.G.cells.num;
            model.thermalConductivity = model.thermalConductivity*ones(nc, 1);
            model.heatCapacity = model.heatCapacity*ones(nc, 1);
            
        end

        function state = updateHeatFlux(model, state)

            k = model.EffectiveThermalConductivity;
            T = state.T;
            
            jHeat = assembleFlux(model, T, k); 

            state.jHeat = jHeat;
            
        end
            
        function state = updateEnergyConservation(model, state)
        % Here, we assume that the fields are updated
        % - jHeatBcSource
        % - jHeatSource 
        % - accumHeat
            
            state = model.updateHeatFlux(state);
            
            flux   = state.jHeat;
            bcsource = state.jHeatBcSource;
            source = state.jHeatSource;
            accum  = state.accumHeat;
            
            energyCons = assembleConservationEquation(model, flux, bcsource, source, accum);
            
            state.energyCons = energyCons;
            
        end
            
        function state = updateOhmSource(model, state)
            
        % The default is a the effective electric conductivity given by the model (but this may be overridden, see for example
        % the electrolyte)
            
            conductivity = model.EffectiveElectronicConductivity;
            state = updateOhmSourceFunc(model, state, conductivity);
            
        end
                    
        function state = updateHeatSource(model, state)
            
            state.jHeatSource = state.jHeatOhmSource; 
            
        end

    end
end

