classdef ThermoElectronicComponent < ElectronicComponent
    
    properties
        
       thermalConductivity
       heatCapacity
       ohmicResistance
       
    end
    
    methods
        
        function model = ThermoElectronicComponent(paramobj)
            
            model = model@ElectronicComponent(paramobj);
            
            fdnames = {'thermalConductivity', ...
                       'heatCapacity', ...
                       'ohmicResistance'};
            model = dispatchParams(model, paramobj, fdnames);
            
            nc = model.G.cells.num;
            model.thermalConductivity = model.thermalConductivity*ones(nc, 1);
            model.heatCapacity = model.heatCapacity*ones(nc, 1);
            model.ohmicResistance = model.ohmicResistance*ones(nc, 1);
        
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
            bcsource = state.jHeatBcSource;
            source = state.jHeatSource;
            accum  = state.accumHeat;
            
            energyCons = assembleConservationEquation(model, flux, bcsource, source, accum);
            
            state.energyCons = energyCons;
            
        end
            
        function state = updateOhmSource(model, state)
            
            state = updateOhmSourceFunc(model, state);
            
        end
                    
        function state = updateHeatSource(model, state)
            
            state.jHeatSource = state.jHeatOhmSource; 
            
        end
                            
    end
end

