classdef OpenElectrode < ElectronicComponent
    
    properties
        
        compInd   % mapping structure for component indices
        phaseInd  % mapping structure for phase indices
        liquidInd % mapping structure for component indices
        gasInd    % mapping structure for component indices

        solidVolumefraction
        leverettCoefficient
        theta % water contact angle
        Permeability

        tau % 
        
        sp % species struct 
        sp.OH.MW
        sp.OH.V0
        sp.OH.D % diffustion coefficient
        sp.OH.t 
        sp.OH.z 
        sp.K.MW
        sp.K.V0
        sp.H2O.MW
        sp.H2O.beta % interpolation coefficient for water equilibrium
        sp.V0 % indexed values for partial molar volumes
        
    end
    
    methods
        
        function model = OpenElectrode(paramobj)

        % paramobj is instance of ElectronicComponentInputParams
            model = model@ElectronicComponent(paramobj);
            
            compInd.H2Oliquid = 1;
            compInd.H2Ogas    = 2;
            compInd.OH        = 3;
            compInd.K         = 4;
            compInd.ncomp     = 5;
            compInd.liquid    = [1; 3; 4];
            compInd.gas       = [2; 5];
            compInd.phaseMap  = [1; 2; 1; 1; 2]; % first component (H2Oliquid) is in phase indexed by 1 (liquid phase), and so on
            
            model.compInd = compInd;
            
            phaseInd.liquid = 1;
            phaseInd.gas    = 2;
            phaseInd.solid  = 3;
            phaseInd.nphase = 3;
            
            model.phaseInd = phaseInd;
            
            liquidInd.H2Oliquid = 1;
            liquidInd.OH = 2;
            liquidInd.K  = 3;
            liquidInd.ncomp  = 3;
            liqudInd.compMap = [1; 3; 4];
            
            model.liquidInd = liquidInd;            
            
            gasInd.H2Ogas  = 1;
            gasInd.ncomp   = 2;
            gasInd.compMap = [2; 5];
            
            model.gasInd = gasInd;
            
            
        end

        function state = updateVolumeFractions(model, state)
            liqeps = state.liqeps;

            state.volumeFractions{model.phaseInd.liquid} = liqeps;
            state.volumeFractions{model.phaseInd.solid} = model.solidVolumefraction;
            state.volumeFractions{model.phaseInd.gas} = 1 - (liqeps + model.solidVolumefraction);
        end


        function state = updateLiquidDensity(model, state)

            state.liquidDensity = state.rholiqeps./state.volumeFractions{model.phaseInd.liquid};
        end
        
        
        function state = updateLiquidPressure(model, state)
        % assemble liquid pressure using capillary pressure function
            K       = model.Permeability;
            levcoef = model.leverettCoefficient;
            theta   = model.theta;
            
            pgas = state.pressures{model.phaseInd.gas};

            vl = state.volumeFractions{model.phaseInd.liquid};
            vg = state.volumeFractions{model.phaseInd.gas};
            vs = state.volumeFractions{model.phaseInd.solid};
            
            % Liquid saturation
            s = vl./(vl + vg);
            
            pc = 0.0694 .* cosd(theta) ./ sqrt(K./vs) .* leverett(levcoef, s);
                        
            pliq = pgas - pc;

            state.pressures{model.phaseInd.liquid} = pliq;
        end
        

        function state = updateOHConcentration(model, state)
            
            vfliquid = state.volumeFractions{model.phaseInd.liquid};
            OHceps = state.OHceps;
            
            state.concentrations{model.liquidInd.OH} = OHceps./vfliquid
        end
        
        
        function state = updateMassH2Ogas(model, state)
        % assemble mass of H2O in gas phase

            state.compGasMasses{compind.H2Ogas} = state.H2Ogasrhoeps;
        end

        function state = updateConcentrations(model, state)
            
            sp = model.sp;

            cOH = state.concentrations{model.liquidInd.OH};
            massliquid = state.liquidDensity;

            cK   = cOH;
            cH2O = (massliquid - cOH.*sp.OH.MW - cK.*sp.K.MW)./sp.H2O.MW;
            cH   = 1e3.*(10.^-sp.H2O.beta .* (1e-3.*cOH).^-1);            
            
            lInd = model.liquidInd;
            state.concentrations{lInd.K}   = cK;
            state.concentrations{lInd.H2O} = cH2O;
            state.concentrations{lInd.H}   = cH;
            
            warning('should update water activity');
            
        end
        

        function state = updateLiquidViscosity(model, state)

            c = state.concentrations(model.liquidInd.OH);
            T = state.T;
            
            mu_par = [4.3e-1  ;
                      -2.51e-2;
                      10^-4   ; 
                      1.3e-1];

            state.viscosities{model.phaseInd.liquid} = exp(mu_par(1) + mu_par(2).*(T - 273.15) + mu_par(3).*(T - 273.15).^2 + mu_par(4).*(1e-3.*c));
            
        end
        
        function state = updateVaporPressure(model, state)

            m = state.OHmolality;
            T = state.T;
            
            a = -1.508e-2    .* m ...
                -1.6788e-3   .* m.^2 ...
                + 2.25887e-5 .* m.^3;
            
            
            b = 1 ...
                - 1.2062e-3 .* m ...
                + 5.6024e-4 .* m.^2;
            
            pw0 = 10.^(35.4462 - (3343.93./T) - 10.9.*log10(T) + 4.1645e-3.*T);
                    
            state.vaporPressure = 10.^(a + b.*log10(pw0)).*1e5;
        end
        
        function state = updateEvaporationTerm(model, state)
            
            %% assemble evaporation term
            psat  = model.sp.H2O;
            MWH2O = model.sp.H2O.MW
            kLV   = 1;
            R     = model.constants.R;
            
            pH2Ovap = state.vaporPressure;
            T       = state.T;
            mH2Ogas = state.compGasMasses{model.gasInd.H2Ogas};
            vl      = state.volumeFractions{model.phaseInd.liquid};
            vg      = state.volumeFractions{model.phaseInd.gas};
            
            pH2Ogas = (mH2Ogas./MWH2O).*R.*T./vg;
            sLiq = vl./(vl + vg);

            evapSrc = 0*p; % initialize r (useful when AD)
            ind = (pH2Ogas > pH2Ovap); 
            if any(ind)
                evapSrc(ind) = (pH2Ovap(ind) - pH2Ogas(ind)) .* kLV ./MW;
            end
            if any(~ind)
                evapSrc(~ind) = sLiq(~ind) .* (pH2Ovap(~ind) - pH2Ogas(~ind)) .* kLV ./ MW;
            end
            
            state.evapH2OSources = evapSrc;
        end
        
        function state = updateViscosities(model, state)
            error('Virtual function. Viscosities are dependent of the particular electrode');
        end
        
        
        function state = updatePhaseVelocities(model, state)
            %% Assemble phase velocities
            K = model.Permeability;
            
            phaseinds = [model.phaseInd.liquid; model.phaseInd.gas];
            for ind = 1 : numel(phaseinds)
                phind = phaseinds{ind};
                p = state.pressures{phind};
                mu = state.viscosities{phind};
                v{ind} = assembleFlux(model, p, K./mu);
            end
            
            state.phaseVelocities = v;
        end
        
        function state = updateOHConvectionFlux(model, state)
            cOH = state.concentrations{model.liquidInd.OH};
            vl  = state.phaseVelocities{model.phaseInd.liquid};
            vfl = state.volumeFractions{model.phaseInd.liquid};
            
            state.convOHFlux = cOH.*vfl.^(1.5).vfl;
        end
        
        
        function state = updateOHDiffusionFlux(model, state)
        % assemble OH diffusion flux
            D = model.sp.OH.D;

            cOH = state.concentrations{model.liqInd.OH};
            vf = state.volumeFractions{model.phaseInd.liquid};
            
            state.diffOHFlux = assembleFlux(model, cOH, vf.^1.5.*D);
        end
        
        function state = updateOHMigrationFlux(model, state)
        % assemble OH migration flux
            F = model.constants.F
            t = model.sp.OH.t;
            z = model.sp.OH.z;
            
            j = state.j;
            
            state.migOHFlux = t./(z.*F).*j;
        end
        
        function state = updateGasFluxes(model, state)
        % assemble convection fluxes

            vfg = state.volumeFractions{model.phaseInd.gas};
            vg  = state.phaseVelocities{model.phaseInd.gas};

            for ind = 1 : model.gasInd.ncomp
                indcomp = model.gasInd.compMap(ind);
                % Note the power 0.5. We use Bruggeman coefficient 1.5 but the component mass is given per *total* volume
                % (meaning that it already is multipled by the volume fraction vf{ph}).
                state.convFluxes{indcomp} =  state.compMasses{indcomp}.*vfg.^0.5).*vg;
            end
            
        end
        
        function state = updateLiquidFlux(model, state)
            liqrho = state.liquidDensity;
            liqvf   = state.volumeFractions{model.phaseInd.liquid};
            liqv    = state.phaseVelocities{model.phaseInd.liquid};

            % We use Bruggeman coefficient 1.5 
            state.liquidFlux = liqrho.*liqvf.^1.5.*liqv;
        end

        function state = updateMolality(model, state)
            
            MW = model.sp.OH.MW;
            
            rho = state.liquidDensity;
            c   = state.concentrations{model.liquidInd.OH};
            
            state.OHmolality = c./(rho - c.*MW);
        end
        
        
        function state = updateESource(model, state)
        % Assemble charge source
            F = model.constants.F;
            z = model.sp.OH.z
            
            OHsrc = state.compOHSources;
            
            state.eSource = OHsrc.*F.*z:
        end
        
        function state = updateGasMassCons(model, state)
        % Assemble mass conservation equations for components in gas phase

            indH2Ogas = model.gasInd.compMap(model.compInd.H2Ogas);
            state.compGasSources{indH2Ogas} = state.compGasSources{indH2Ogas} + state.evapH2Osources;
            
            for ind = 1 : model.gasInd.ncomp
                indcomp = model.gasInd.compMap(ind);
                gasMassCons{ind} = assembleConservationEquation(model, ...
                                                                state.compGasFluxes{ind}   , ...
                                                                state.compGasBcSources{ind}, ...
                                                                state.compGasSources{ind}  , ...
                                                                state.accumTerms{indcomp});
            end
            
            state.gasMassCons = gasMassCons;
        end

        function state = updateLiquidMassCons(model, state)
        % Assemble mass conservation for the overall liquid component
            
            liquidSource = state.OHsource.*(model.sp.OH.MW + model.sp.K.MW) ... 
                           + state.H2OliquidSource.*model.sp.H2O.MW ...
                           - state.evapH2Osource;

            state.liquiMassCons = assembleConservationEquation(model          , ...
                                                              state.liquidFlux, ...
                                                              liquidSource    , ...
                                                              0               , ...
                                                              state.liquidAccumTerm);
        end
        
        function state = updateOHMassCons(model, state)
        % Assemble mass conservation equation for OH

            OHflux = state.convOHFlux + state.diffOHFlux + state.migOHFlux;
            state.OHMassCons = assembleConservationEquation(model         , ...
                                                            OHFlux        , ...
                                                            state.OHsource, ...
                                                            state.accumTerms{model.compInd.OH});
        end

        function state = updatePartialMolarVolumes(model, state)
            error('not used in first implementation');
            OH = model.sp.OH;
            K = model.sp.K;
            
            m   = state.OHmolality;
            rho = state.rhoLiquid;
            T   = state.T;
            
            rho_par = [778.6106  ;
                       42.8554   ;  
                       1.7400    ;
                       -1.4428   ;
                       0.04585   ;
                       -0.0034   ;
                       0.0197    ;
                       4.1087e-4 ;
                       -1.1139e-4];
                        
            % Empirical model for first derivative of mass density with respect to KOH molality
            drhodm = rho_par(2) + ...
                     2.*rho_par(4).*m + rho_par(5).*T + ...
                     3.*rho_par(7).*m.^2 + ...
                     2.*rho_par(8).*m.*T + rho_par(9).*T.^2;
            
            pmv = -1./(rho.^2).*drhodm.*(1 + m.*MW) + (1./rho).*MW;
            
            state.partialMolarVolumes{model.liquidInd.OH} = pmv.*abs(OH.V0)./(abs(OH.V0) + abs(K.V0));
            state.partialMolarVolumes{model.liquidInd.K} = pmv.*abs(K.V0)./(abs(OH.V0) + abs(K.V0));
        end
        
        
        function state = liquidStateEquation(model, state)
            
            liqStateEq = -1;
            for ind = 1 : model.liquidInd.ncomp
                liqStateEq = liqStateEq + state.concentrations{ind}.*model.sp.V0(ind);
            end
            
            state.liquidStateEquation = liqStateEq;
        end
        
    end


end

