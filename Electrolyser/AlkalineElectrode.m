classdef AlkalineElectrode < ElectronicComponent
    
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
        % sp.H2O.mu0 % Standard chemical potential

        
        % sp.V0 % indexed values for partial molar volumes
        
    end
    
    methods
        
        function model = AlkalineElectrode(paramobj)

            % paramobj is instance of ElectronicComponentInputParams
            paramobj.use_thermal = false;
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
            % total liquid density (mass of liquid per total volume)
            varnames{end + 1} = 'liqrhoeps';
            % Phase pressures
            phasePressures = VarName({}, 'phasePressures', nph);
            varnames{end + 1} = phasePressures;
            % Phase volume fractions
            volumeFractions = VarName({}, 'volumeFractions', nph);
            varnames{end + 1} = volumeFractions;
            % Masses for each component of gas (per total volume, as used in mass of conservation law)
            varnames{end + 1}  = VarName({}, 'compGasPressures', ngas);
            % Masses for each component of gas (per total volume, as used in mass of conservation law)
            compGasMasses = VarName({}, 'compGasMasses', ngas);
            varnames{end + 1} = compGasMasses;
            % Liquid density (Mass of liquid per volume of liquid)
            varnames{end + 1} = 'liqrho';
            % Concentrations in the liquid
            concentrations = VarName({}, 'concentrations', nliquid);
            varnames{end + 1} = concentrations;
            % Concentrations in the OH molality
            % varnames{end + 1} = 'OHmolality';
            % H2O activity (we could assemble all activities there but it seems that only H2O activity is needed)
            varnames{end + 1} = 'H2Oa';
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
            
            % Vapor pressure
            varnames{end + 1} = 'vaporPressure';
            
            
            %% Coupling variables
            
            % Mass sources for the gas Components
            varnames{end + 1}  = VarName({}, 'compGasSources', ngas);
            % Mass sources at the boundaries for the gas components
            varnames{end + 1}  = VarName({}, 'compGasBcSources', ngas);
            % Accumulation term for the gass components in [kg/s]
            varnames{end + 1}  = VarName({}, 'compGasAccums', ngas);
            % Source of OH (in mole)
            varnames{end + 1} = 'OHSource';
            % Liquid-Vapor exchange rate for H2O (H2Oliquid <-> H2Ogas)
            varnames{end + 1} = 'H2OliquidVaporExchangeRate';
            % Source of H2Oliquid (in mole)
            varnames{end + 1} = 'H2OliquidSource';
            % Accumulation terms for OH in [mol/s]
            varnames{end + 1} = 'OHaccum';
            % Accumulation term for the overall liquid components in [kg/s]
            varnames{end + 1} = 'liquidAccumTerm';
            
            %% Residual variables            
            % Residual for the mass conservation equations of the components in the gas
            varnames{end + 1}  = VarName({}, 'compGasMassCons', ngas);
            % Residual for the mass conservation equation of the aggregated components in the liquid phase
            varnames{end + 1} = 'liquidMassCons';
            % Residual for the conservation equation for OH in the liquid phase (unit is Mol)
            varnames{end + 1} = 'OHMassCons';
            % Residual for the equation of state of the liquid
            varnames{end + 1} = 'liquidStateEquation';

            model = model.registerVarNames(varnames);


            fn = @() AlkalineElectrode.updateVolumeFractions;
            inputnames = {'liqeps'};
            model = model.registerPropFunction({volumeFractions, fn, inputnames});            

            % assemble liquid pressure using capillary pressure function
            fn = @() AlkalineElectrode.updateLiquidPressure;
            inputnames = {VarName({}, 'phasePressures', nph, phaseInd.gas), volumeFractions};
            model = model.registerPropFunction({VarName({}, 'phasePressures', nph, phaseInd.liquid), fn, inputnames});
            
            %% assemble masses 
            
            % assemble mass of OH
            fn = @() AlkalineElectrode.updateOHconcentration;
            inputnames = {'OHceps', ...
                          VarName({}, 'volumeFractions', nph, phaseInd.liquid)};
            var = VarName({}, 'concentrations', nliquid, liquidInd.OH);
            model = model.registerPropFunction({var, fn, inputnames});            
            
            % assemble mass of H2O in gas phase
            fn = @() AlkalineElectrode.updateMassH2Ogas;
            inputnames = {'H2Ogasrhoeps'};
            var = VarName({}, 'compGasMasses', ngas, gasInd.H2Ogas);
            model = model.registerPropFunction({var, fn, inputnames});            
            
            % update liquid density
            fn = @() AlkalineElectrode.updateLiquidDensity;
            inputnames = {'liqrhoeps', ...
                          VarName({}, 'volumeFractions', nph, phaseInd.liquid)};
            model = model.registerPropFunction({'liqrho', fn, inputnames});            
            
            % assemble concentrations
            fn = @() AlkalineElectrode.updateConcentrations;
            inputnames = {'liqrho', ...
                          VarName({}, 'concentrations', nliquid, liquidInd.OH)};
            ind = setdiff([1 : nliquid]', liquidInd.OH);
            model = model.registerPropFunction({VarName({}, 'concentrations', nliquid, ind), fn, inputnames});            
            
            fn = @() AlkalineElectrode.updateWaterActivity;
            inputnames = {VarName({}, 'compGasPressures', ngas, gasInd.H2Ogas)};
            model = model.registerPropFunction({'H2Oa', fn, inputnames});

            % compute OH molalities
            fn = @() AlkalineElectrode.updateMolality;
            inputnames = {VarName({}, 'concentrations', nph, model.liquidInd.OH), 'liqrho'};
            model = model.registerPropFunction({'OHmolality', fn, inputnames});            
            
            %% assemble vapor pressure
            fn = @() AlkalineElectrode.updateVaporPressure;
            inputnames = {'T', 'OHmolality'};
            model = model.registerPropFunction({'vaporPressure', fn, inputnames});
            
            %% assemble evaporation term
            fn = @() AlkalineElectrode.updateEvaporationTerm;
            inputnames = {'T', ...
                          'vaporPressure', ...
                          VarName({}, 'compGasMasses', ngas, gasInd.H2Ogas), ...
                          volumeFractions};
            model = model.registerPropFunction({'H2OliquidVaporExchangeRate', fn, inputnames});

            %% Assemble Liquid viscosity
            % We move this to specific electrode
            % see [[file:~/Matlab/Projects/AlkalineElectrolyzerContinuumModel/Model/Materials/O2mix.m::function viscosity(obj)][for the gas]] 
            % see [[file:~/Matlab/Projects/AlkalineElectrolyzerContinuumModel/Model/Materials/H2mix.m::function viscosity(obj)][for the gas]]
            % see [[file:~/Matlab/Projects/AlkalineElectrolyzerContinuumModel/Model/Materials/KOH.m::function viscosity(obj)][for the liquid]]
            
            % fn = @() AlkalineElectrode.updateLiquidViscosity;
            % inputnames = {'T', ...
                          % VarName({}, 'concentrations', nliquid, liquidInd.OH)};
            % model = model.registerPropFunction({VarName({}, 'viscosities', nph, phaseInd.liquid ), fn, inputnames});
            
            if ~model.useZeroDmodel
                
                %% Assemble phase velocities
                fn = @() AlkalineElectrode.updatePhaseVelocities;
                
                inputnames = {VarName({}, 'phasePressures', nph, phaseInd.mobile), ...
                              VarName({}, 'viscosities', nph, phaseInd.mobile)};
                model = model.registerPropFunction({VarName({}, 'phaseVelocities', nph, phaseInd.mobile), fn, inputnames});
                
                
                %% Assemble OH specific fluxes
                
                % assemble OH convection flux
                fn = @() AlkalineElectrode.updateOHConvectionFlux;
                inputnames = {VarName({}, 'concentrations', nliquid, liquidInd.OH), ...
                              VarName({}, 'phaseVelocities', nph, phaseInd.liquid), ...
                              VarName({}, 'volumeFractions', nph, phaseInd.liquid)};
                model = model.registerPropFunction({'convOHFlux', fn, inputnames}); 
                
                % assemble OH diffusion flux
                fn = @() AlkalineElectrode.updateOHDiffusionFlux;
                inputnames = {VarName({}, 'concentrations', nliquid, liquidInd.OH), ...
                              VarName({}, 'volumeFractions', nph, phaseInd.liquid)};
                model = model.registerPropFunction({'diffOHFlux', fn, inputnames}); 
                
                % assemble OH migration flux
                fn = @() AlkalineElectrode.updateOHMigrationFlux;
                inputnames = {'j'};
                model = model.registerPropFunction({'migOHFlux', fn, inputnames}); 
                
                %% Assemble  fluxes
                
                % assemble fluxes of gas components
                fn = @() AlkalineElectrode.updateGasFluxes;
                inputnames = {compGasMasses, ...
                              VarName({}, 'volumeFractions', nph, phaseInd.gas), ...
                              VarName({}, 'phaseVelocities', nph, phaseInd.gas)};
                model = model.registerPropFunction({VarName({}, 'compGasFluxes', ngas), fn, inputnames});
                
                % Assemble flux of the overall liquid component
                fn = @() AlkalineElectrode.updateLiquidFlux;
                inputnames = {VarName({}, 'phaseVelocities', nph, phaseInd.liquid), ...
                              'liqrho', ...
                              VarName({}, 'volumeFractions', nph, phaseInd.liquid)};
                model = model.registerPropFunction({'liquidFlux', fn, inputnames});
                
            end
            
            %% Assemble charge source

            fn = @() AlkalineElectrode.updateESource;
            inputnames = {'OHSource'};
            model = model.registerPropFunction({'eSource', fn, inputnames});
            
            %% Assemble the residual equations


            if model.useZeroDmodel
                % Assemble mass conservation equations for components in gas phase
                fn = @() AlkalineElectrode.updateGasMassCons0;

                for igas = 1 : ngas
                    inputnames = {VarName({}, 'compGasBcSources', ngas, igas), ...
                                  VarName({}, 'compGasSources'  , ngas, igas), ...
                                  VarName({}, 'compGasAccums'   , ngas, igas)};
                    if igas == gasInd.H2Ogas
                        inputnames{end + 1} = 'H2OliquidVaporExchangeRate';
                    end
                    model = model.registerPropFunction({VarName({}, 'compGasMassCons', nph, igas), fn, inputnames});
                end


                % Assemble mass conservation for the overall liquid component
                fn = @() AlkalineElectrode.updateLiquidMassCons0;
                inputnames = {'liquidAccumTerm', ...
                              'OHSource'       , ...
                              'H2OliquidSource', ...
                              'H2OliquidVaporExchangeRate'};
                model = model.registerPropFunction({'liquidMassCons', fn, inputnames});
                
                % Assemble mass conservation equation for OH
                fn = @() AlkalineElectrode.updateOHMassCons0;
                inputnames = {'OHSource', ...
                              'OHaccum'};                          
                model = model.registerPropFunction({'OHMassCons', fn, inputnames});

            else

                fn = @() AlkalineElectrode.updateLiquidViscosity;
                inputnames = {'T', ....
                              VarName({}, 'concentrations', nliquid, liquidInd.OH)};
                model = model.registerPropFunction({VarName({}, 'viscosities', nph, phaseInd.liquid), fn, inputnames});
                
                % Assemble mass conservation equations for components in gas phase
                fn = @() AlkalineElectrode.updateGasMassCons;
                for igas = 1 : ngas
                    inputnames = {VarName({}, 'compGasBcSources', ngas, igas), ...
                                  VarName({}, 'compGasFluxes'   , ngas, igas), ...
                                  VarName({}, 'compGasSources'  , ngas, igas), ...
                                  VarName({}, 'compGasAccums'   , ngas, igas)};
                    if igas == gasInd.H2Ogas
                        inputnames{end + 1} = 'H2OliquidVaporExchangeRate';
                    end
                    model = model.registerPropFunction({VarName({}, 'compGasMassCons', nph, igas), fn, inputnames});
                end

                % Assemble mass conservation for the overall liquid component
                fn = @() AlkalineElectrode.updateLiquidMassCons;
                inputnames = {'liquidFlux'     , ...
                              'liquidAccumTerm', ...
                              'OHSource'       , ...
                              'H2OliquidSource', ...
                              'H2OliquidVaporExchangeRate'};
                model = model.registerPropFunction({'liquidMassCons', fn, inputnames});
                
                % Assemble mass conservation equation for OH
                fn = @() AlkalineElectrode.updateOHMassCons;
                inputnames = {'convOHFlux', ...
                              'diffOHFlux', ...
                              'migOHFlux' , ...
                              'OHSource'  , ...
                              'OHaccum'};
                model = model.registerPropFunction({'OHMassCons', fn, inputnames});

            end
            
            % Assemble partial Molar Volumes (not used in first implementation)
            % fn = @() AlkalineElectrode.updatePartialMolarVolumes;
            % inputnames = {'OHmolality', 'T', 'liqrho'};
            % ind = [model.liquidInd.OH; model.liquidInd.K];
            % model = model.registerPropFunction({VarName({}, 'partialMolarVolumes', nliquid, ind), fn, inputnames});

            % we use liquid incompressibility for the moment
            fn = @() AlkalineElectrode.updateLiquidAccum;
            inputnames = {};
            model = model.registerPropFunction({'liquidAccumTerm', fn, inputnames});
            
            % Assemble residual of equation of state for the liquid phase
            fn = @() AlkalineElectrode.liquidStateEquation;
            inputnames = {concentrations};
            model = model.registerPropFunction({'liquidStateEquation', fn, inputnames});

            fn = @() AlkalineElectrode.updateAccumTerms;
            functionCallSetupFn = @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction);
            fn = {fn, functionCallSetupFn};
            inputnames = {'H2Ogasrhoeps'};
            model = model.registerPropFunction({VarName({}, 'compGasAccums', ngas, gasInd.H2Ogas), fn, inputnames});
            inputnames = {'OHceps'};
            model = model.registerPropFunction({'OHaccum', fn, inputnames});

            model = model.removeVarName(VarName({}, 'phasePressures', nph, phaseInd.solid));
            
        end
        

        function state = updateVolumeFractions(model, state)

            liqeps = state.liqeps;

            state.volumeFractions{model.phaseInd.liquid} = liqeps;
            state.volumeFractions{model.phaseInd.solid}  = model.solidVolumefraction;
            state.volumeFractions{model.phaseInd.gas}    = 1 - (liqeps + model.solidVolumefraction);
            
        end


        function state = updateLiquidDensity(model, state)

            vf = state.volumeFractions{model.phaseInd.liquid};

            state.liqrho = state.liqrhoeps./vf;
            
        end
        
        
        function state = updateLiquidPressure(model, state)
        % assemble liquid pressure using capillary pressure function

            K       = model.Permeability;
            levcoef = model.leverettCoefficient;
            theta   = model.theta;
            
            pgas = state.phasePressures{model.phaseInd.gas};

            vl = state.volumeFractions{model.phaseInd.liquid};
            vg = state.volumeFractions{model.phaseInd.gas};
            vs = state.volumeFractions{model.phaseInd.solid};
            
            % Liquid saturation
            s = vl./(vl + vg);
            
            pc = 0.0694 .* cosd(theta) ./ sqrt(K./vs) .* leverett(levcoef, s);
                        
            pliq = pgas - pc;

            state.phasePressures{model.phaseInd.liquid} = pliq;
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
            massliquid = state.liqrho;

            cK   = cOH;
            cH2O = (massliquid - cOH.*sp.OH.MW - cK.*sp.K.MW)./sp.H2O.MW;
            cH   = 1e3.*(10.^-sp.H2O.beta .* (1e-3.*cOH).^-1);            
            
            lInd = model.liquidInd;
            state.concentrations{lInd.K}   = cK;
            state.concentrations{lInd.H2O} = cH2O;
            state.concentrations{lInd.H}   = cH;
            
        end

        function state = updateWaterActivity(model, state)

            R = model.constants.R;
            mu0 = model.sp.mu0;
            
            p = state.compGasPressures(model.gasInd.H2O);
            T = state.T;
            
            state.H2Oa = mu0 + R*T*log(p ./ 1e5);
            
        end
        
        function state = updateLiquidViscosity(model, state)

            cOH = state.concentrations{model.liquidInd.OH};
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
            
            state.H2OliquidVaporExchangeRate = evapSrc;
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
                p = state.phasePressures{phind};
                mu = state.viscosities{phind};
                v{ind} = assembleFlux(model, p, K./mu);
            end
            
            state.phaseVelocities = v;
        end
        
        function state = updateOHConvectionFlux(model, state)
            
            cOH = state.concentrations{model.liquidInd.OH};
            v  = state.phaseVelocities{model.phaseInd.liquid};
            vf = state.volumeFractions{model.phaseInd.liquid};

            state.convOHFlux = cOH.*vf.^(1.5).*v;
            
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

            vf = state.volumeFractions{model.phaseInd.gas};
            v  = state.phaseVelocities{model.phaseInd.gas};

            for igas = 1 : model.gasInd.ncomp
                % NOTE: We take the power to 0.5 of the volume fraction but it corresponds to a Bruggeman coefficient
                % 1.5, because the component mass is given per *total* volume
                % (meaning that it already is multipled by the volume fraction vf).
                state.compGasFluxes{igas} =  state.compGasMasses{igas}.*(vf.^0.5).*v;
            end
            
        end
        
        function state = updateLiquidFlux(model, state)
             
            liqrho = state.liqrho;
            liqvf  = state.volumeFractions{model.phaseInd.liquid};
            liqv   = state.phaseVelocities{model.phaseInd.liquid};

            % We use Bruggeman coefficient 1.5 
            state.liquidFlux = liqrho.*liqvf.^1.5.*liqv;
            
        end

        function state = updateMolality(model, state)
            
            MW = model.sp.OH.MW;
            
            rho = state.liqrho;
            c   = state.concentrations{model.liquidInd.OH};
            
            state.OHmolality = c./(rho - c.*MW);
            
        end
        
        
        function state = updateESource(model, state)
        % Assemble charge source

            F = model.constants.F;
            z = model.sp.OH.z
            
            OHsrc = state.OHSource;
            
            state.eSource = OHsrc.*F.*z;
            
        end
        

        function state = updateGasMassCons(model, state)
        % Assemble mass conservation equations for components in gas phase

            indH2Ogas = model.gasInd.compMap(model.compInd.H2Ogas);
            state.compGasSources{indH2Ogas} = state.compGasSources{indH2Ogas} + state.evapH2Osources;
            
            for ind = 1 : model.gasInd.ncomp
                compGasMassCons{ind} = assembleConservationEquation(model, ...
                                                                state.compGasFluxes{ind}   , ...
                                                                state.compGasBcSources{ind}, ...
                                                                state.compGasSources{ind}  , ...
                                                                state.compGasAccums{ind});
            end
            
            state.compGasMassCons = compGasMassCons;
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
                                                            state.OHaccum);
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

