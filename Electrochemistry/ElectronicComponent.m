classdef ElectronicComponent < PhysicalModel
    
    properties
        
        constants; % physical constants, see PhysicalConstants
        EffectiveElectronicConductivity;
        
    end

    methods
        
        function model = ElectronicComponent(paramobj)
        % paramobj is instance of ElectronicComponentInputParams
            model = model@PhysicalModel([]);

            % OBS : All the models should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            fdnames = {'G', ...
                       'EffectiveElectronicConductivity'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
            % setup discrete differential operators
            model.operators = localSetupOperators(model.G);
            
            model.constants = PhysicalConstants();
            
        end

        function state = updateCurrent(model, state)
            
            sigmaeff = model.EffectiveElectronicConductivity;
            phi = state.phi;
            
            j = assembleFlux(model, phi, sigmaeff); 

            state.j = j;
            
        end
        
        function state = updateChargeConservation(model, state)
            
            state = model.updateCurrent(state);

            flux   = state.j;
            bcsource = state.jBcSource;
            source = state.eSource;
            accum  = 0;
            
            chargeCons = assembleConservationEquation(model, flux, bcsource, source, accum);
            
            state.chargeCons = chargeCons;
            
        end
        
    end
end

