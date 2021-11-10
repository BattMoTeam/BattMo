classdef Electrolyte < ElectroChemicalComponent
    
    properties
        
        sp
        compnames
        ncomp

        indchargecarrier
        chargeCarrierName
        
        Separator
        
        volumeFraction
        
        thermalConductivity % intrinsic thermal conductivity value
        heatCapacity % intrinsic heat capacity value

        EffectiveThermalConductivity
        EffectiveHeatCapacity

        updateConductivityFunc
        updateDiffusionCoefficientFunc
        
    end

    methods
        
        function model = Electrolyte(paramobj)
        % paramobj is instance of ElectrolyteInputParams or a derived class
            model = model@ElectroChemicalComponent(paramobj);
            
            model.Separator = Separator(paramobj.Separator);
            
            fdnames = {'sp'                 , ...
                       'compnames'          , ...
                       'chargeCarrierName'  , ...
                       'thermalConductivity', ...
                       'heatCapacity'};

            model = dispatchParams(model, paramobj, fdnames);

            model.updateConductivityFunc = str2func(paramobj.updateConductivityFunc.functionname);
            model.updateDiffusionCoefficientFunc = str2func(paramobj.updateDiffusionCoefficientFunc.functionname);            

            model.ncomp = numel(model.compnames);
            [isok, indchargecarrier] = ismember(model.chargeCarrierName, model.compnames);
            assert(isok, 'charge carrier not found in the list of components');
            model.indchargecarrier = indchargecarrier;
            
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
        
        function state = updateConcentrations(model, state)
            
            cs{1} = state.c;
            cs{2} = state.c;
            
            state.cs = cs;
            
        end
                

        function state = updateConductivity(model, state)
            
            cLi = state.cs{1};
            T = state.T;
            
            func = model.updateConductivityFunc;
            
            state.conductivity = func(cLi, T);
        end

        function state = updateChemicalCurrent(model, state)
            
            cLi          = state.cs{1}; % concentration of Li+
            T            = state.T;     % temperature
            phi          = state.phi;   % potential
            conductivity = state.conductivity;   % potential
            
            cs  = state.cs;         
            
            ncomp = model.ncomp; % number of components
            sp = model.sp;

            % calculate the concentration derivative of the chemical potential for each species in the electrolyte
            R = model.constants.R;
            for ind = 1 : ncomp
                dmudcs{ind} = R .* T ./ cs{ind};
            end
                        
            % volume fraction of electrolyte
            volfrac = model.volumeFraction;
            % Compute effective ionic conductivity in porous media
            conductivityeff = conductivity .* volfrac .^1.5;
            
            % setup chemical fluxes
            jchems = cell(1, ncomp);
            F = model.constants.F;
            for i = 1 : ncomp
                coeff = conductivityeff .* sp.t(i) .* dmudcs{i} ./ (sp.z(i).*F);
                jchems{i} = assembleFlux(model, cs{i}, coeff);
            end
            
            state.dmudcs = dmudcs;
            state.jchems = jchems;
            
        end
        
        function state = updateDiffusionCoefficient(model, state)
            
            c = state.cs{1};
            T = state.T;

            func = model.updateDiffusionCoefficientFunc;
            
            D = func(c, T);
            
            % set effective coefficient
            state.D = D .* model.volumeFraction .^1.5;
        
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
            
            % volume fraction of electrolyte
            volfrac = model.volumeFraction;
            % Compute effective ionic conductivity in porous media
            conductivityeff = conductivity.*volfrac.^1.5;
            
            j = assembleFlux(model, phi, conductivityeff);
            for ind = 1 : ncomp
                j = j + jchems{ind};
            end

            state.j = j;
            
        end
        
        function state = updateMassFlux(model, state)
            
            ind = model.indchargecarrier;
            
            % We assume that the current and the diffusion coefficient D has been updated when this function is called
            c = state.cs{ind};
            j = state.j;
            D = state.D;
            
            %% 1. Flux from diffusion
            diffFlux = assembleFlux(model, c, D);
            state.diffFlux = diffFlux;
            
            %% 2. Flux from electrical forces
            F = model.constants.F;
            fluxE = model.sp.t(ind) ./ (model.sp.z(ind) .* F) .* j;
            
            %% 3. Sum the two flux contributions
            flux = diffFlux + fluxE;
            
            state.massFlux = flux;
            
        end
        
    end
end

