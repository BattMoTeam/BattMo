classdef ThermoElectronicComponent < ElectronicComponent
    
    properties
        
       thermalConductivity
       heatCapacity
       ohmicResistance
       
    end
    
    methods
        
        function model = ThermoElectronicComponent(paramobj)
            
            model = model@ElectronicComponent(component);
            
            fdnames = {'thermaConductivity', ...
                       'heatCapacity', ...
                       'ohmicResistance'};
            model = dispatchParams(model, paramobj, fdnames);
            
        end

        function state = updateHeatFlux(model, state)

            k = model.thermalConductivity;
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
            bcflux = state.jHeatBcSource;
            source = state.jHeatSource;
            accum  = state.accumHeat;
            
            energyCons = assembleConservationEquation(model, flux, bcflux, source, accum);
            
            state.chargeCons = chargeCons;
            
        end
            
        function state = updateOhmSource(model, state)
            
            state = updateOhmSourceFunc(model, state);
            
        end
                    
        function state = updateHeatSource(model, state)
            
            state.jHeatSource = state.jHeatOhmSource; 
            
        end
                            
    end
end

