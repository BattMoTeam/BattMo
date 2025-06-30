classdef LithiumPlatingLatz < BaseModel
% Implementation of the model from Hein et al
% @article{hein2020electrochemical,
%   title={An electrochemical model of lithium plating and stripping in lithium ion batteries},
%   author={Hein, Simon and Danner, Timo and Latz, Arnulf},
%   journal={ACS Applied Energy Materials},
%   volume={3},
%   number={9},
%   pages={8519--8531},
%   year={2020},
%   publisher={ACS Publications}
% }
    
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
        platedConcentrationRef                  % Reference concentration for plated lithium [mol/m^3]
        c
        volumetricSurfaceArea                   % Interfacial surface area between graphite and electrolyte per unit volume [m²/m³]
        particleRadius                          % Radius of graphite particles, used in solid-state diffusion modeling
        volumeFraction                          % Porosity
        
        platedLiMonolayerThickness              % Thickness of one monolayer of plated lithium around a particle [m]. See equation (S-3)

        SEIFraction                             % Fraction of graphite surface covered by SEI (solid electrolyte interphase)
        MSEI                                    % Molar mass of SEI [kg/mol]
        rhoSEI                                  % Density of SEI [kg/m³]
        deltaSEI0                               % Initial thickness of SEI layer [m]
        sigmaSEI                                % Ionic conductivity of SEI [S/m]

        useSEI                                  % Boolean flag: whether SEI effects are included in overpotential calculations

        % numerical parameters
        logReg  = 1e-6 % Value used for regularization of the logarithm in the
                       % activity of plated lithium, to avoid numerical issues when
                       % the concentration is very low (default 1e-6)
        
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
                       'platedConcentrationRef', ...
                       'volumetricSurfaceArea' , ...
                       'particleRadius',...
                       'volumeFraction', ...
                       'platedLiMonolayerThickness', ...
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

            varnames{end + 1} = 'T';                        % Temperature [K]

            varnames{end + 1} = 'OCP';                      % Open circuit Voltage [V]

            varnames{end + 1} = 'phiElectrode';             % Potential of the solid electrode [V]

            varnames{end + 1} = 'phiElectrolyte';           % Potential of the electrolyte [V]

            varnames{end + 1} = 'cElectrolyte';             % Concentration in Li of the electrolyte [mol/m^3]

            varnames{end + 1} = 'platedConcentrationNorm';  % Normalized plated Li concentration [unitless]

            varnames{end + 1} = 'platedConcentration';      % Actual plated concentration [mol/m^3]

            varnames{end + 1} = 'platedConcentrationAccum'; % Accumulated plated concentration [mol/m^3]

            varnames{end + 1} = 'platedConcentrationCons';  % Conservation equation for plated concentration [mol/m^2/s]

            % update this comment
            
            % Warning : we use the conventions from Rein et al paper, where the fluxes are considered from the inside to
            % the outside, except for the intercalation flux which has been defined before for which the convention is
            % that the flux goes from the electrolyte to the electrode

            varnames{end + 1} = 'platingFlux';       % Plating flux [mol/m^2/s], platingFlux is positive if stripping and negative if plating

            varnames{end + 1} = 'chemicalFlux';      % Flux of plated lithium into the electrode [mol/m^2/s], chemicalFlux is Negative if platedLithium inserts in the electrode

            varnames{end + 1} = 'etaPlating';        % Overpotential for plated lithium (Butler-Volmer) [V]

            varnames{end + 1} = 'etaChemical';       % Overpotential for plated insertion (Butler-Volmer) [V]

            varnames{end + 1} = 'activityPlated';    % Activity of plated lithium [unitless]

            varnames{end + 1} = 'surfaceCoverage';   % Fraction of Surface covered by plated lithium [unitless]

            varnames{end + 1} = 'cElectrodeSurface'; % Concentration of solid lithium at Surface of particule [mol/m3]

            % varnames{end + 1} = 'platedThickness';   % Thickness of plated lithium around a particle [m]. See equation (27)

            model = model.registerVarNames(varnames);

            fn = @LithiumPlatingLatz.updatePlatedConcentrationAccum;
            fn = {fn, @(pf) PropFunction.accumFuncCallSetupFn(pf)};
            model = model.registerPropFunction({'platedConcentrationAccum', fn, {'platedConcentration'}});

            fn = @LithiumPlatingLatz.updateActivityPlated;
            model = model.registerPropFunction({'activityPlated', fn, {'platedConcentration'}});

            fn = @LithiumPlatingLatz.updatePlatedConcentration;
            model = model.registerPropFunction({'platedConcentration', fn, {'platedConcentrationNorm'}});

            fn = @LithiumPlatingLatz.updateEtaPlating;
            model = model.registerPropFunction({'etaPlating', fn, {'phiElectrode', 'phiElectrolyte', 'activityPlated', 'T'}});
            
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
            
            % fn = @LithiumPlatingLatz.updatePlatedThickness;
            % model = model.registerPropFunction({'platedThickness', fn, {'platedConcentration'}});

        end
        
        function state = updatePlatedConcentration(model, state)
            state.platedConcentration = state.platedConcentrationNorm * model.platedConcentrationRef;
        end

        function state = updatePlatedConcentrationAccum(model, state, state0, dt)
            state.platedConcentrationAccum = (state.platedConcentration - state0.platedConcentrationNorm * model.platedConcentrationRef) / dt;
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

        function state = updateEtaPlating(model, state)
            
            phiS = state.phiElectrode;
            phiE = state.phiElectrolyte;
            aPl  = state.activityPlated;
            T    = state.T;

            eta = phiS - phiE  + (model.R * T / model.F) .* log(aPl + model.logReg);
            state.etaPlating = eta;
        end

        function state = updateEtaChemical(model, state)
            
            aPl = state.activityPlated;
            T   = state.T;
            OCP = state.OCP;

            eta = -OCP -(model.R * T / model.F) .* log(aPl + model.logReg);
            state.etaChemical = eta;
            
        end

        function state = updatePlatingFlux(model, state)

            R = model.R;
            F = model.F;
            
            eta = state.etaPlating;
            ce  = state.cElectrolyte;
            T   = state.T;

            th = model.platedConcentrationRef * 0.1;

            i0 = model.kPl * regularizedPow(ce, model.alphaPl, th);
            j = i0 .* (exp((model.alphaPl * F * eta) / (R * T)) - ...
                       exp((-model.alphaStr * F * eta) / (R * T)));
            
            state.platingFlux = j; %mol/s/m2
            
        end

        function state = updateChemicalFlux(model, state)

            eta = state.etaChemical;
            T   = state.T;
            F = model.F;
            cSo = state.cElectrodeSurface;

            th = model.platedConcentrationRef * 0.001;
            jCh = model.kChInt * regularizedSqrt(cSo, th)  * (exp(model.alphaChInt * F * eta / (model.R * T)) - ...
                                  exp(-model.alphaChInt * F * eta / (model.R * T)));
            
            state.chemicalFlux = jCh; %mol/s/m2
            
        end

        function state = updatePlatedConcentrationCons(model, state)

            vsa = model.volumetricSurfaceArea;
            F   = model.F;
            
            surfcov = state.surfaceCoverage;
            
            % a part of the flux goes in the SEI if the option is active
            % So only 1 - model.SEIFraction goes in the particles
            % coeff_plating = 1 - model.SEIFraction * model.useSEI;
            % flux  = (state.chemicalFlux - state.platingFlux*coeff_plating) * s;

            flux  = (state.platingFlux - state.chemicalFlux) * surfcov;
            
            accum = state.platedConcentrationAccum;
            
            % multiplying the flux by vsa to get a volumetric flux
            % (mol/s/m3)
            
            state.platedConcentrationCons = accum + flux*vsa;
            
        end

        function state = updateSurfaceCoverage(model, state)
            
            nLimit = model.nPlLimit; % n of plated lithium necessary to cover the whole surface of the particle
            r      = model.particleRadius;
            vf     = model.volumeFraction;
            
            cLimit = nLimit * vf / ((4/3)*pi*r^3); %and we computed the associated concentration

            platedConcentration = state.platedConcentration;
            % check for AD
            state.surfaceCoverage = min(platedConcentration ./ cLimit, 1.0);
            
        end
        % 
        % function state = updatePlatedThickness(model, state)
        % 
        %     nLimit = model.nPlLimit; % n of plated lithium necessary to cover the whole surface of the particle
        %     r      = model.particleRadius;
        %     vf     = model.volumeFraction;
        %     thck1layer = model.platedLiMonolayerThickness;
        % 
        %     cLimit = nLimit * vf / ((4/3)*pi*r^3); %and we computed the associated concentration
        %     cPl = state.platedConcentration;
        % 
        %     dPl = thck1layer .* cPl ./ cLimit;
        % 
        %     state.platedThickness = dPl;
        % 
        % end
    end
end

