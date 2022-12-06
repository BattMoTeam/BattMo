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

        useZeroDmodel
        
        sp % species struct 
        % sp.OH.MW
        % sp.OH.V0
        % sp.OH.D % diffustion coefficient
        % sp.OH.t 
        % sp.OH.z 
        % sp.K.MW
        % sp.K.V0
        % sp.H2O.MW
        % sp.H2O.beta % interpolation coefficient for water equilibrium
        % sp.V0 % indexed values for partial molar volumes
        
    end
    
    methods
        
        function model = OpenElectrode(paramobj)

        % paramobj is instance of ElectronicComponentInputParams
            model = model@ElectronicComponent(paramobj);
            
            compInd.H2Oliquid = 1;
            compInd.H2Ogas    = 2;
            compInd.OH        = 3;
            compInd.K         = 4;
            % compInd.(H2 or O2) = 5 % should be instantiated by derived class see HydrogenElectrode.m and OxygenElectrode.m
            compInd.ncomp     = 5;
            compInd.liquid    = [1; 3; 4];
            compInd.gas       = [2; 5];
            compInd.phaseMap  = [1; 2; 1; 1; 2]; % first component (H2Oliquid) is in phase indexed by 1 (liquid phase), and so on
            
            model.compInd = compInd;
            
            phaseInd.liquid = 1;
            phaseInd.gas    = 2;
            phaseInd.solid  = 3;
            phaseInd.mobile = [1; 2];
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

        function model = registerVarAndPropfuncNames(model)
            
        %% Declaration of the Dynamical Variables and Function of the model
        % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@ElectronicComponent(model);

            phaseInd  = model.phaseInd;
            liquidInd = model.liquidInd;
            gasInd    = model.gasInd;
            compInd   = model.compInd;

            ncomp   = compInd.ncomp;
            ngas    = gasInd.ncomp;
            nph     = phaseInd.nphase;
            nliquid = liquidInd.ncomp;
            nmobph  = numel(phaseInd.mobile);

            varnames = {};
            % Total concentration of OH- (mol per total volume, that is not only liquid volume)
            varnames{end + 1} = 'OHceps';
            % Total density of gas H2O  (mass per total volume, that is not only gas volume)
            varnames{end + 1} = 'H2Ogasrhoeps';
            % Liquid volume fraction
            varnames{end + 1} = 'liqeps';
            % Phase pressures
            pressures = VarName({}, 'pressures', nph);
            varnames{end + 1} = pressures;
            % Phase volume fractions
            volumeFractions = VarName({}, 'volumeFractions', nph);
            varnames{end + 1} = volumeFractions;
            
            % Masses for each component of gas (per total volume, as used in mass of conservation law)
            % compGasMasses = VarName({}, 'compGassMasses', ngas);
            % varnames{end + 1} = compGasMasses;

            % Liquid density (Mass of liquid per volume of liquid)
            varnames{end + 1} = 'liquidDensity';
            % Concentrations in the liquid
            concentrations = VarName({}, 'concentrations', nliquid);
            varnames{end + 1} = concentrations;
            % Concentrations in the OH molality
            % varnames{end + 1} = 'OHmolality';
            % H2O activity (we could assemble all activities there but it seems that only H2O activity is needed)
            varnames{end + 1} = 'H2Oactivity';
            % OH molar volume
            % partialMolarVolumes = VarName({}, 'partialMolarVolumes', nliquid);
            % varnames{end + 1} = partialMolarVolumes;


            if ~model.useZeroDmodel
                
                %% Flow variables

                viscosities = VarName({}, 'viscosities', nmobph);
                varnames{end + 1} = viscosities;            
                
                %% Phase velocities

                phaseVelocities = VarName({}, 'phaseVelocities', nmobph);
                varnames{end + 1} = phaseVelocities;
                
                %% Fluxes

                % Mass fluxes for the gass components
                compGasFluxes = VarName({}, 'compGasFluxes', ngas);
                varnames{end + 1} = compGasFluxes;
                % Convective fluxes for OH  
                varnames{end + 1} = 'convOHFlux';
                % Diffusion flux for OH
                varnames{end + 1} = 'diffOHFlux';
                % Migration flux for OH
                varnames{end + 1} = 'migOHFlux';
                % Mass flux for total of liquid components
                varnames{end + 1} = 'liquidFlux';            

            end
            
            %% Phase transfer variables
            
            % Vapor pressure
            varnames{end + 1} = 'vaporPressure';
            % Evaporation sources for the water (H2Oliquid <-> H2Ogas)
            varnames{end + 1} = 'evapH2OSource';
            
            %% Coupling variables
            
            % Mass sources for the gas Components
            compGasSources = VarName({}, 'compGasSources', ngas);
            varnames{end + 1} = compGasSources;
            % Mass sources at the boundaries for the gas components
            compGasBcSources = VarName({}, 'compGasBcSources', ngas);
            varnames{end + 1} = compGasBcSources;
            % Source of OH (in mole)
            varnames{end + 1} = 'OHSource';
            % Source of H2Oliquid (in mole)
            varnames{end + 1} = 'H2OliquidSource';
            % Accumulation terms (for all components)
            % beware of the units (typically in kg but in mol for OH)
            accumTerms = VarName({}, 'accumTerms', ncomp);
            varnames{end + 1} = accumTerms;
            % Accumulation term for the overall liquid components (mass)
            varnames{end + 1} = 'liquidAccumTerm';
            
            %% Residual variables            
            % Residual for the mass conservation equations of the components in the gas
            gasMassCons = VarName({}, 'gasMassCons', ngas);
            varnames{end + 1} = gasMassCons;
            % Residual for the mass conservation equation of the aggregated components in the liquid phase
            varnames{end + 1} = 'liquidMassCons';
            % Residual for the conservation equation for OH in the liquid phase (unit is Mol)
            varnames{end + 1} = 'OHMassCons';
            % Residual for the equation of state of the liquid
            varnames{end + 1} = 'liquidStateEquation';

            model = model.registerVarNames(varnames);


            fn = @() OpenElectrode.updateVolumeFractions;
            inputnames = {'liqeps'};
            model = model.registerPropFunction({volumeFractions, fn, inputnames});            

            % assemble liquid pressure using capillary pressure function
            fn = @() OpenElectrode.updateLiquidPressure;
            inputnames = {VarName({}, 'pressures', nph, phaseInd.gas), volumeFractions};
            model = model.registerPropFunction({VarName({}, 'pressures', nph, phaseInd.liquid), fn, inputnames});
            
            %% assemble masses 
            
            % assemble mass of OH
            fn = @() OpenElectrode.updateOHconcentration;
            inputnames = {'OHceps', ...
                          VarName({}, 'volumeFractions', nph, phaseInd.liquid)};
            var = VarName({}, 'concentrations', nliquid, liquidInd.OH);
            model = model.registerPropFunction({var, fn, inputnames});            
            
            % assemble mass of H2O in gas phase
            fn = @() OpenElectrode.updateMassH2Ogas;
            inputnames = {'H2Ogasrhoeps'};
            var = VarName({}, 'compGasMasses', ngas, gasInd.H2Ogas);
            model = model.registerPropFunction({var, fn, inputnames});            
            
            % update liquid density
            fn = @() OpenElectrode.updateLiquidDensity;
            inputnames = {'rholiqeps', ...
                          VarName({}, 'volumeFractions', nph, phaseInd.liquid)};
            model = model.registerPropFunction({'liquidDensity', fn, inputnames});            
            
            % assemble concentrations
            fn = @() OpenElectrode.updateConcentrations;
            inputnames = {'liquidDensity', ...
                          VarName({}, 'concentrations', nliquid, liquidInd.OH)};
            ind = setdiff([1 : nliquid]', liquidInd.OH);
            model = model.registerPropFunction({VarName({}, 'concentrations', nliquid, ind), fn, inputnames});            
            model = model.registerPropFunction({'H2Oactivity', fn, inputnames});

            % compute OH molalities
            fn = @() OpenElectrode.updateMolality;
            inputnames = {VarName({}, 'concentrations', nph, model.liquidInd.OH), 'liquidDensity'};
            model = model.registerPropFunction({'OHmolality', fn, inputnames});            
            
            %% assemble vapor pressure
            fn = @() OpenElectrode.updateVaporPressure;
            inputnames = {'T', 'OHmolality'};
            model = model.registerPropFunction({'vaporPressure', fn, inputnames});
            
            %% assemble evaporation term
            fn = @() OpenElectrode.updateEvaporationTerm;
            inputnames = {'T', ...
                          'vaporPressure', ...
                          VarName({}, 'compGasMasses', ngas, gasInd.H2Ogas), ...
                          volumeFractions};
            model = model.registerPropFunction({'evapH2OSource', fn, inputnames});

            %% Assemble Liquid viscosity
            % We move this to specific electrode
            % see [[file:~/Matlab/Projects/AlkalineElectrolyzerContinuumModel/Model/Materials/O2mix.m::function viscosity(obj)][for the gas]] 
            % see [[file:~/Matlab/Projects/AlkalineElectrolyzerContinuumModel/Model/Materials/H2mix.m::function viscosity(obj)][for the gas]]
            % see [[file:~/Matlab/Projects/AlkalineElectrolyzerContinuumModel/Model/Materials/KOH.m::function viscosity(obj)][for the liquid]]
            
            % fn = @() OpenElectrode.updateLiquidViscosity;
            % inputnames = {'T', ...
                          % VarName({}, 'concentrations', nliquid, liquidInd.OH)};
            % model = model.registerPropFunction({VarName({}, 'viscosities', nph, phaseInd.liquid ), fn, inputnames});
            
            if ~model.useZeroDmodel
                
                %% Assemble phase velocities
                fn = @() OpenElectrode.updatePhaseVelocities;
                
                inputnames = {VarName({}, 'pressures', nph, phaseInd.mobile), ...
                              VarName({}, 'viscosities', nph, phaseInd.mobile)};
                model = model.registerPropFunction({VarName({}, 'phaseVelocities', nph, phaseInd.mobile), fn, inputnames});
                
                
                %% Assemble OH specific fluxes
                
                % assemble OH convection flux
                fn = @() OpenElectrode.updateOHConvectionFlux;
                inputnames = {VarName({}, 'concentrations', nliquid, liquidInd.OH), ...
                              VarName({}, 'phaseVelocities', nph, phaseInd.liquid), ...
                              VarName({}, 'volumeFractions', nph, phaseInd.liquid)};
                model = model.registerPropFunction({'convOHFlux', fn, inputnames}); 
                
                % assemble OH diffusion flux
                fn = @() OpenElectrode.updateOHDiffusionFlux;
                inputnames = {VarName({}, 'concentrations', nliquid, liquidInd.OH), ...
                              VarName({}, 'volumeFractions', nph, phaseInd.liquid)};
                model = model.registerPropFunction({'diffOHFlux', fn, inputnames}); 
                
                % assemble OH migration flux
                fn = @() OpenElectrode.updateOHMigrationFlux;
                inputnames = {'j'};
                model = model.registerPropFunction({'migOHFlux', fn, inputnames}); 
                
                %% Assemble  fluxes
                
                % assemble fluxes of gas components
                fn = @() OpenElectrode.updateGasFluxes;
                inputnames = {'compGasMasses', ...
                              VarName({}, 'volumeFractions', nph, phaseInd.gas), ...
                              VarName({}, 'phaseVelocities', nph, phaseInd.gas)};
                getGasInd = @(ind) VarName({}, 'compGasFluxes', ncomp, ind);
                model = model.registerPropFunction({getGasInd(1), fn, inputnames});
                model = model.registerPropFunction({getGasInd(2), fn, inputnames});
                
                % Assemble flux of the overall liquid component
                fn = @() OpenElectrode.updateLiquidFlux;
                inputnames = {VarName({}, 'phaseVelocities', nph, phaseInd.liquid), ...
                              'liquidDensity', ...
                              VarName({}, 'volumeFractions', nph, phaseInd.liquid)};
                model = model.registerPropFunction({'liquidFlux', fn, inputnames});
                
            end
            
            %% Assemble charge source

            fn = @() OpenElectrode.updateESource;
            inputnames = {'OHSource'};
            model = model.registerPropFunction({'eSource', fn, inputnames});
            
            %% Assemble the residual equations


            if model.useZeroDmodel
                % Assemble mass conservation equations for components in gas phase
                fn = @() OpenElectrode.updateGasMassCons0;
                inputnames = {'evapH2OSource' , ...
                              VarName({}, 'compGasSources', ngas, gasInd.H2Ogas), ...
                              VarName({}, 'accumTerms', ncomp, compInd.H2Ogas)};
                model = model.registerPropFunction({VarName({}, 'gasMassCons', nph, gasInd.H2Ogas), fn, inputnames});

                % Assemble mass conservation for the overall liquid component
                fn = @() OpenElectrode.updateLiquidMassCons0;
                inputnames = {'liquidAccumTerm', ...
                              'OHSource'       , ...
                              'H2OliquidSource', ...
                              'evapH2OSource'};
                model = model.registerPropFunction({'liquidMassCons', fn, inputnames});
                
                % Assemble mass conservation equation for OH
                fn = @() OpenElectrode.updateOHMassCons0;
                inputnames = {'OHSource'  , ...
                              VarName({}, 'accumTerms', ncomp, compInd.OH)};                          
                model = model.registerPropFunction({'OHMassCons', fn, inputnames});

            else
                
                warning('update from useZeroDmodel above')
                % Assemble mass conservation equations for components in gas phase
                fn = @() OpenElectrode.updateGasMassCons;
                inputnames = {'evapH2OSource' , ...
                              compGasSources  , ...
                              compGasFluxes   , ...
                              compGasBcSources, ...
                              VarName({}, 'accumTerms', ncomp, compInd.gas)};
                model = model.registerPropFunction({'gasMassCons', fn, inputnames});

                % Assemble mass conservation for the overall liquid component
                fn = @() OpenElectrode.updateLiquidMassCons;
                inputnames = {'liquidFlux'     , ...
                              'liquidAccumTerm', ...
                              'OHSource'       , ...
                              'H2OliquidSource', ...
                              'evapH2OSource'};
                model = model.registerPropFunction({'liquidMassCons', fn, inputnames});
                
                % Assemble mass conservation equation for OH
                fn = @() OpenElectrode.updateOHMassCons;
                inputnames = {'convOHFlux', ...
                              'diffOHFlux', ...
                              'migOHFlux' , ...
                              'OHSource'  , ...
                              VarName({}, 'accumTerms', ncomp, compInd.OH)};                          
                model = model.registerPropFunction({'OHMassCons', fn, inputnames});

            end
            
            % Assemble partial Molar Volumes (not used in first implementation)
            % fn = @() OpenElectrode.updatePartialMolarVolumes;
            % inputnames = {'OHmolality', 'T', 'liquidDensity'};
            % ind = [model.liquidInd.OH; model.liquidInd.K];
            % model = model.registerPropFunction({VarName({}, 'partialMolarVolumes', nliquid, ind), fn, inputnames});

            % we use liquid incompressibility for the moment
            fn = @() OpenElectrode.updateLiquidAccum;
            inputnames = {};
            model = model.registerPropFunction({'liquidAccumTerm', fn, inputnames});
            
            % Assemble residual of equation of state for the liquid phase
            fn = @() OpenElectrode.liquidStateEquation;
            inputnames = {concentrations};
            model = model.registerPropFunction({'liquidStateEquation', fn, inputnames});

            fn = @() OpenElectrode.updateAccumTerms;
            inputnames = {'H2Ogasrhoeps'};
            model = model.registerPropFunction({VarName({}, 'accumTerms', ncomp, compInd.H2Ogas), fn, inputnames});
            inputnames = {'OHceps'};
            model = model.registerPropFunction({VarName({}, 'accumTerms', ncomp, compInd.OH), fn, inputnames});
            inputnames = {'OHceps'};
            model = model.registerPropFunction({VarName({}, 'accumTerms', ncomp, compInd.OH), fn, inputnames});

            model = model.registerStaticVarName('T');
            
            model = model.removeVarName(VarName({}, 'pressures', nph, phaseInd.solid));
            model = model.removeVarName(VarName({}, 'accumTerms', ncomp, [compInd.H2Oliquid, compInd.K]));
            
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

            warning('move this function to specific electrode as updateGasViscosity')
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

            state.convOHFlux = cOH.*vfl.^(1.5).*vfl;
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
                state.convFluxes{indcomp} =  state.compMasses{indcomp}.*(vfg.^0.5).*vg;
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
            
            state.eSource = OHsrc.*F.*z;
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

        function state = updateGasMassCons0(model, state)
        % zero flux version
            for ind = 1 : model.gasInd.ncomp
                state.compGasFluxes{ind} = 0;
            end
            
            state = model.updateGasMassCons(state);
            
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

        function state = updateLiquidMassCons0(model, state)
        % zero flux version
            state.liquidFlux = 0;
            state = model.updateLiquidMassCons(state)
        end

        
        function state = updateOHMassCons(model, state)
        % Assemble mass conservation equation for OH

            OHflux = state.convOHFlux + state.diffOHFlux + state.migOHFlux;
            state.OHMassCons = assembleConservationEquation(model         , ...
                                                            OHFlux        , ...
                                                            state.OHsource, ...
                                                            state.accumTerms{model.compInd.OH});
        end

        function state = updateOHMassCons0(model, state)
        % zero flux version

            state.convOHFlux = 0;
            state.diffOHFlux = 0;
            state.migOHFlux  = 0;
            state = model.updateOHMassCons(state);

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

