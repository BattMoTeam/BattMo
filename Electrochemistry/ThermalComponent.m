classdef ThermalComponent < PhysicalModel
    
    properties
        
       EffectiveThermalConductivity
       EffectiveHeatCapacity % in [J][K]^-1[m]^-3

       couplingTerm
       externalHeatTransferCoefficient
       externalTemperature
       
    end
    
    methods
        
        function model = ThermalComponent(paramobj)
            
            model = model@PhysicalModel([]);
            
            % OBS : All the models should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            fdnames = {'G', ...
                       'EffectiveThermalConductivity', ...
                       'EffectiveHeatCapacity', ...
                       'couplingTerm', ...
                       'externalHeatTransferCoefficient', ...
                       'externalTemperature'};
            model = dispatchParams(model, paramobj, fdnames);
            
            % setup discrete differential operators
            model.operators = localSetupOperators(model.G);
            
            if ~isempty(model.EffectiveThermalConductivity)
                nc = model.G.cells.num;
                model.EffectiveThermalConductivity = model.EffectiveThermalConductivity*ones(nc, 1);
                model.EffectiveHeatCapacity = model.EffectiveHeatCapacity*ones(nc, 1);
            end
        
        end

        function state = updateHeatFlux(model, state)

            k = model.EffectiveThermalConductivity;
            T = state.T;
            
            jHeat = assembleFlux(model, T, k); 

            state.jHeat = jHeat;
            
        end
        
        function state = updateHeatSourceTerm(model, state)
        % sum up the heat source terms
            state.jHeatSource = state.jHeatOhmSource + state.jHeatChemicalSource + state.jHeatReactionSource;
            
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

        function state = updateThermalBoundarySourceTerms(model, state)

            G          = model.G;
            coupterm   = model.couplingTerm;
            T_ext      = model.externalTemperature;
            lambda_ext = model.externalHeatTransferCoefficient;
            lambda     = model.EffectiveThermalConductivity;
            
            coupcells = coupterm.couplingcells;
            coupfaces = coupterm.couplingfaces;
            nc = model.G.cells.num;
              
            T = state.T;
            T = T(coupcells);
            
            if isempty(coupfaces)
                % 1D case (External faces are not available, we consider a volumetric cooling instead)
                A = G.cells.volumes(coupcells);
                t_eff = lambda_ext*A;
            else
                % Face couling (multidimensional case)
                A = G.faces.areas(coupfaces);
                t_ext = lambda_ext.*A;
                
                t = model.operators.harmFaceBC(lambda, coupfaces);
                
                t_eff = 1./(1./t + 1./t_ext);
            end
            
            jHeatBcSource = zeros(nc, 1);
            if isa(T, 'ADI')
                adsample = getSampleAD(T);
                adbackend = model.AutoDiffBackend;
                jHeatBcSource = adbackend.convertToAD(jHeatBcSource, adsample);
            end
            
            jHeatBcSource(coupcells) = t_eff.*(T_ext - T);
            
            state.jHeatBcSource = jHeatBcSource;
            
        end
        
    end
    
end

