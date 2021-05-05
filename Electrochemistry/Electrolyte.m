classdef Electrolyte < ElectroChemicalComponent
    
    properties
        
        indchargecarrier
        
        sp
        compnames
        ncomp

        Separator
        
        volumeFraction
        
                            
        thermalConductivity % intrinsic thermal conductivity value
        heatCapacity % intrinsic heat capacity value

        EffectiveThermalConductivity
        EffectiveHeatCapacity
        
    end

    methods
        
        function model = Electrolyte(paramobj)
        % paramobj is instance of ElectrolyteInputParams or a derived class
            model = model@ElectroChemicalComponent(paramobj);
            
            model.Separator = Separator(paramobj.sep);
            
            fdnames = {'sp'                 , ...
                       'compnames'          , ...
                       'indchargecarrier'   , ...
                       'ncomp'              , ...
                       'thermalConductivity', ...
                       'heatCapacity'};

            model = dispatchParams(model, paramobj, fdnames);
            
            sep = 'Separator';

            % We set the electrolyte volumeFraction based on the porisity of the separator
            G = model.G;
            Gp = G.mappings.parentGrid;
            
            model.volumeFraction = NaN(G.cells.num, 1);
            
            elyte_cells = zeros(Gp.cells.num, 1);
            elyte_cells(G.mappings.cellmap) = (1 : model.G.cells.num)';
            elyte_cells_sep = elyte_cells(model.(sep).G.mappings.cellmap);
            model.volumeFraction(elyte_cells_sep) = model.(sep).porosity;
            
            % The effective thermal conductivity in the common region between electrode and electrolyte is setup when the battery is
            % set up. Here we set up the thermal conductivity of the electrolyte in the separator region (we assume for
            % the moment constant values for both porosity and thermal conductivity but this can be changed).
            
            model.EffectiveThermalConductivity = NaN(G.cells.num, 1);
            model.EffectiveThermalConductivity(elyte_cells_sep) = model.(sep).porosity.*model.thermalConductivity;
            
            
        end
        
        function state = updateChemicalCurrent(model, state)
            error('virtual function');
        end
        
        function state = updateDiffusionCoefficient(model, state)
            error('virtual function');
        end
        
        function state = updateCurrentBcSource(model, state)
        % no boundary current fluxes (only volumetric from the reactions)
            state.jBcSource = 0;
        end
        
        function state  = updateCurrent(model, state) 
           
            ncomp = model.ncomp;
            
            phi = state.phi;
            jchems = state.jchems;
            conductivity = state.conductivity;
            
            j = assembleFlux(model, phi, conductivity);
            for ind = 1 : ncomp
                j = j + jchems{ind};
            end

            state.j = j;
            
        end
        
        function state = updateChargeCarrierFlux(model, state)
            
            ind = model.indchargecarrier;
            
            % We assume that the current and the diffusion coefficient D has been updated when this function is called
            c = state.cs{ind};
            j = state.j;
            D = state.D;
            
            %% 1. Flux from diffusion
            fluxDiff = assembleFlux(model, c, D);
            
            %% 2. Flux from electrical forces
            F = model.constants.F;
            fluxE = model.sp.t{ind} ./ (model.sp.z{ind} .* F) .* j;
            
            %% 3. Sum the two flux contributions
            flux = fluxDiff + fluxE;
            
            %% 4. Apply scaling (maybe not the right place but consistent  with assembleConservationEquation - at
            %% least when this comment has beem written...)
            flux = flux*F; 
            
            state.LiFlux = flux;
            
        end
        
    end
end

