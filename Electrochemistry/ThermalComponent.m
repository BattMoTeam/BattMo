classdef ThermalComponent < PhysicalModel
    
    properties
        
       EffectiveThermalConductivity
       EffectiveHeatCapacity % in [J][K]^-1[m]^-3
       
    end
    
    methods
        
        function model = ThermalComponent(paramobj)
            
            model = model@PhysicalModel([]);
            
            % OBS : All the models should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            fdnames = {'EffectiveThermalConductivity', ...
                       'EffectiveHeatCapacity'};
            model = dispatchParams(model, paramobj, fdnames);
            
            % setup discrete differential operators
            model.operators = localSetupOperators(model.G);
            
            model.constants = PhysicalConstants();
            
            nc = model.G.cells.num;
            model.EffectiveThermalConductivity = model.EffectiveThermalConductivity*ones(nc, 1);
            model.EffectiveHeatCapacity = model.EffectiveHeatCapacity*ones(nc, 1);
        
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
            
            flux     = state.jHeat;
            bcsource = state.jHeatBcSource;
            source   = state.jHeatSource;
            accum    = state.accumHeat;
            
            energyCons = assembleConservationEquation(model, flux, bcsource, source, accum);
            
            state.energyCons = energyCons;
            
        end
            
    end
    
end

