classdef LithiumPlatingLatz < BaseModel

    properties

        F = PhysicalConstants.F
        R = PhysicalConstants.R

        alphaPl                                 % Symmetry factor for lithium plating/stripping reaction (typically ≈ 0.3)
        alphaStr                                % Symmetry factor for stripping (may be 1-alphaPl)
        alphaChInt                              % Symmetry factor for charge-neutral chemical intercalation of plated lithium

        kPl                                     % Reaction rate constant for lithium plating (N⁰₀ Plating)
        kChInt                                  % Reaction rate constant for chemical intercalation of plated lithium
        kInter                                  % Reaction rate constant for direct lithium-ion intercalation into graphite

        nPl0                                    % Phenomenological parameter: minimum lithium amount needed to activate metal activity (see eqn (8))
        nPlLimit                                % Limit amount of plated lithium corresponding to one monolayer on graphite surface (see eqn (26))

        volumetricSurfaceArea                   % Interfacial surface area between graphite and electrolyte per unit volume [m²/m³]
        particleRadius                          % Radius of graphite particles, used in solid-state diffusion modeling
        volumeFraction                          % Porosity

        SEIFraction                             % Fraction of graphite surface covered by SEI (solid electrolyte interphase)
        MSEI                                    % Molar mass of SEI [kg/mol]
        rhoSEI                                  % Density of SEI [kg/m³]
        deltaSEI0                               % Initial thickness of SEI layer [m]
        sigmaSEI                                % Ionic conductivity of SEI [S/m]

        useSEI                                  % Boolean flag: whether SEI effects are included in overpotential calculations

    end

    methods

        function model = LithiumPlatingLatz(inputparams)
            model = model@BaseModel();
            fdnames = {'alphaPl'    , ...    
                       'alphaStr'   , ...   
                       'alphaChInt' , ... 
                       'kPl'        , ...    
                       'kChInt'     , ...
                       'kInter'     , ...
                       'nPl0'       , ...    
                       'nPlLimit'   , ... 
                       'volumetricSurfaceArea' , ...
                       'particleRadius',...
                       'volumeFraction', ...
                       'SEIFraction', ... 
                       'MSEI'       , ...        
                       'rhoSEI'     , ...      
                       'deltaSEI0'  , ...   
                       'sigmaSEI'   , ...    
                       'useSEI'};
            model = dispatchParams(model, inputparams, fdnames);

        end

        function model = registerVarAndPropfuncNames(model)

            useSEI = model.useSEI;

            varnames = {};

            varnames{end + 1} = 'T';                      % Temperature [K]

            varnames{end + 1} = 'OCP';                    % Open circuit Voltage [V]

            varnames{end + 1} = 'phiElectrode';           % Potential of the solid electrode [V]

            varnames{end + 1} = 'phiElectrolyte';         % Potential of the electrolyte [V]

            varnames{end + 1} = 'cElectrolyte';           % Concentration in Li of the electrolyte [mol/m^3]

            varnames{end + 1} = 'platedConcentration';    % Plated Li concentration in the electrode [mol/m^3]

            varnames{end + 1} = 'platedConcentrationAccum'; % Accumulated plated concentration [mol/m^3]

            varnames{end + 1} = 'platedConcentrationCons'; % Conservation equation for plated concentration [mol/m^2/s]

            % Warning : fluxes are considered from the inside to the outside

            % platingFlux is positive if stripping and negative if plating
            % chemicalFlux is Negative if platedLithium inserts in the electrode

            varnames{end + 1} = 'platingFlux';           % Plating flux [C/m^2/s]

            varnames{end + 1} = 'chemicalFlux';          % Flux of plated lithium into the electrode [C/m^2/s]

            varnames{end + 1} = 'etaPlating';            % Overpotential for plated lithium (Butler-Volmer) [V]

            varnames{end + 1} = 'etaChemical';           % Overpotential for plated insertion (Butler-Volmer) [V]

            varnames{end + 1} = 'activityPlated';        % Activity of plated lithium [unitless]

            varnames{end + 1} = 'surfaceCoverage';       % Fraction of Surface covered by plated lithium [unitless]

            varnames{end + 1} = 'cElectrodeSurface';       % Concentration of solid lithium at Surface of particule [mol/m3]

            if useSEI
                
                varnames{end + 1} = 'nSEIAccum';         % Accumulated SEI amount [mol/m^2]

                varnames{end + 1} = 'nSEI';              % SEI amount [mol/m^2]

                varnames{end + 1} = 'nSEICons';          % SEI Conservation term [mol/m^2/s]

                varnames{end + 1} = 'SEIThickness';      % Thickness of the SEI layer [m]
            end

            model = model.registerVarNames(varnames);

            fn = @LithiumPlatingLatz.updatePlatedConcentrationAccum;
            fn = {fn, @(pf) PropFunction.accumFuncCallSetupFn(pf)};
            model = model.registerPropFunction({'platedConcentrationAccum', fn, {'platedConcentration'}});

            fn = @LithiumPlatingLatz.updateActivityPlated;
            model = model.registerPropFunction({'activityPlated', fn, {'platedConcentration'}});

            if useSEI
                fn = @LithiumPlatingLatz.updateEtaPlatingSEI;
                model = model.registerPropFunction({'etaPlating', fn, {'phiElectrode', 'phiElectrolyte', 'activityPlated', 'SEIThickness', 'T'}});
            else
                fn = @LithiumPlatingLatz.updateEtaPlating;
                model = model.registerPropFunction({'etaPlating', fn, {'phiElectrode', 'phiElectrolyte', 'activityPlated', 'T'}});
            end
            fn = @LithiumPlatingLatz.updateEtaChemical;
            model = model.registerPropFunction({'etaChemical', fn, {'activityPlated', 'T', 'OCP'}});

            fn = @LithiumPlatingLatz.updatePlatingFlux;
            model = model.registerPropFunction({'platingFlux', fn, {'cElectrolyte', 'etaPlating', 'T'}});

            fn = @LithiumPlatingLatz.updateChemicalFlux;
            model = model.registerPropFunction({'chemicalFlux', fn, {'etaChemical', 'T', 'cElectrodeSurface'}});

            fn = @LithiumPlatingLatz.updatePlatedConcentrationCons;
            model = model.registerPropFunction({'platedConcentrationCons', fn, {'platedConcentrationAccum', 'platingFlux', 'chemicalFlux', 'surfaceCoverage'}});

            fn = @LithiumPlatingLatz.updateSurfaceCoverage;
            model = model.registerPropFunction({'surfaceCoverage', fn, {'platedConcentration'}});

            if useSEI
                
                fn = @LithiumPlatingLatz.updateNSEIAccum;
                fn = {fn, @(pf) PropFunction.accumFuncCallSetupFn(pf)};
                model = model.registerPropFunction({'nSEIAccum', fn, {'nSEI'}});

                fn = @LithiumPlatingLatz.updateNSEICons;
                model = model.registerPropFunction({'nSEICons', fn, {'platingFlux', 'nSEIAccum'}});

                fn = @LithiumPlatingLatz.updateSEIThickness;
                model = model.registerPropFunction({'SEIThickness', fn, {'nSEI'}});

            end
        end

        function state = updatePlatedConcentrationAccum(model, state, state0, dt)
            state.platedConcentrationAccum = (state.platedConcentration - state0.platedConcentration) / dt;
        end

        function state = updateActivityPlated(model, state)
            platedConcentration = state.platedConcentration;
            
            n0 = model.nPl0;
            r = model.particleRadius;
            poros = model.volumeFraction;
            % switching the n of lithium plated for one particle to a concentration with this volume
            c0 = n0 * poros / ((4/3)*pi*r^3);
            
            state.activityPlated = platedConcentration^4 ./ (platedConcentration^4 + c0^4);            
        end

        % function state = updateEtaPlatingSEI(model, state)
        % 
        %     phiS = state.phiElectrode;
        %     phiE = state.phiElectrolyte;
        %     aPl = state.activityPlated;
        % 
        %     RSEI = state.SEIThickness / model.sigmaSEI;
        %     j = 5e-12;  % Warning warning warning
        % 
        %     eta = phiS - phiE - j * RSEI + (model.R * state.T / model.F) .* log(aPl);
        %     state.etaPlating = eta;
        % end

        function state = updateEtaPlating(model, state)
            
            phiS = state.phiElectrode;
            phiE = state.phiElectrolyte;
            aPl  = state.activityPlated;
            T    = state.T;

            eta = phiS - phiE  + (model.R * T / model.F) .* log(aPl);
            state.etaPlating = eta;
        end

        function state = updateEtaChemical(model, state)
            
            aPl = state.activityPlated;
            T   = state.T;
            OCP = state.OCP;

            eta = -OCP -(model.R * T / model.F) .* log(aPl);
            state.etaChemical = eta;
            
        end

        function state = updatePlatingFlux(model, state)

            R = model.R;
            F = model.F;
            eta = state.etaPlating;
            ce  = state.cElectrolyte;
            T   = state.T;

            i0 = model.kPl * ce.^model.alphaPl;
            j = i0 .* (exp((model.alphaPl * F * eta) / (R * T)) - ...
                       exp((-model.alphaStr * F * eta) / (R * T)));
            
            state.platingFlux = j * F; %C/s/m2
            
        end

        function state = updateChemicalFlux(model, state)

            eta = state.etaChemical;
            T   = state.T;
            F = model.F;
            cSo = state.cElectrodeSurface;
            jCh = model.kChInt * (cSo)^(1/2)  * (exp(model.alphaChInt * F * eta / (model.R * T)) - ...
                                  exp(-model.alphaChInt * F * eta / (model.R * T)));
            
            state.chemicalFlux = jCh * F; %C/s/m2
            
        end

        function state = updatePlatedConcentrationCons(model, state)

            vsa = model.volumetricSurfaceArea;
            F = model.F;
            s = state.surfaceCoverage;
            
            %a part of the flux goes in the SEI if the option is active
            %So only 1 - model.SEIFraction goes in the particles
            coeff_plating = 1 - model.SEIFraction * model.useSEI;
            
            % converting the current into a mass flux with /F

            flux  = (state.chemicalFlux - state.platingFlux*coeff_plating) * s / F;
            accum = state.platedConcentrationAccum;
            
            % multiplying the flux by vsa to get a volumetric flux
            % (mol/s/m3)
            state.platedConcentrationCons = assembleConservationEquation(model, 0, 0, flux*vsa , accum);
            
        end

        function state = updateSurfaceCoverage(model, state)
            
            nLimit = model.nPlLimit; % n of plated lithium necessary to cover the whole surface of the particle
            r = model.particleRadius;
            poros = model.volumeFraction;
            cLimit = nLimit * poros / ((4/3)*pi*r^3); %and we computed the associated concentration

            platedConcentration = state.platedConcentration;
            state.surfaceCoverage = min(platedConcentration ./ cLimit, 1.0);
            
        end

        function state = updateNSEIAccum(model, state, state0, dt)
            state.nSEIAccum = (state.nSEI - state0.nSEI) / dt;
        end

        function state = updateNSEICons(model, state)
            
            flux  = state.platingFlux / F  * model.SEIFraction;
            accum = state.nSEIAccum;
            
            state.nSEICons = assembleConservationEquation(model, 0, 0, flux, accum);
            
        end

        function state = updateSEIThickness(model, state) %need to modify this
            
            delta = model.deltaSEI0 + (model.MSEI * state.nSEI) / (model.rhoSEI * model.F);
            state.SEIThickness = delta;
        end

    end
end
