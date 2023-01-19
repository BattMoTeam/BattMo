classdef PorousTransportLayer < ElectronicComponent
    
    properties


        Boundary
        
        compInd     % mapping structure for component indices
        phaseInd    % mapping structure for phase indices
        mobPhaseInd % mapping structure for mobile phase indices
        liquidInd   % mapping structure for component indices
        gasInd      % mapping structure for component indices

        solidVolumeFraction  % Solid volume fraction
        leverettCoefficients % coefficients for the Leverett function definition, see method updateLiquidPressure and function leverett.m
        theta                % water contact angle
        permeability         % Permeability [Darcy]
        BruggemanCoefficient % Bruggeman coefficient (not used yet)
        
        sp % species struct
        % sp.OH.MW    : Molecular weight of OH [kg mol^-1]
        % sp.OH.V0    : Partial molar volume of OH [m^3 mol^-1]
        % sp.OH.D     : Diffusion coefficient [m^2 s^-1]
        % sp.OH.t     : Transference coefficient [-]
        % sp.OH.z     : Charge number [-]
        % sp.K.MW     : Molecular weight of K [kg mol^-1]
        % sp.K.V0     : Partial molar volume of K [m^3 mol^-1]
        % sp.H2O.MW   : Molecular weight of H2O [kg mol^-1]
        % sp.H2O.kLV  : Liquid-vapor exchange rate
        % sp.H2O.mu0  : Standard chemical potential
        % sp.H2O.V0   : Partial molar volume of H2O [m^3 mol^-1]        
        
        V0s % indexed values for partial molar volumes (helper that is setup at initialization)

        externalCouplingTerm
        
    end
    
    methods
        
        function model = PorousTransportLayer(paramobj)
            
            % paramobj is instance of ElectronicComponentInputParams
            paramobj.use_thermal = false;
            model = model@ElectronicComponent(paramobj);

            %% Setup the model using the input parameters
            fdnames = {'G'                   , ...
                       'solidVolumeFraction' , ... 
                       'leverettCoefficients', ... 
                       'theta'               , ... 
                       'BruggemanCoefficient', ... 
                       'permeability'        , ...
                       'sp'                  , ...
                       'externalCouplingTerm'};
            model = dispatchParams(model, paramobj, fdnames);

            
            compInd.H2Oliquid = 1;
            compInd.H2Ogas    = 2;
            compInd.OH        = 3;
            compInd.K         = 4;
            compInd.activeGas = 5;
            % compInd.(H2 or O2) = 5 % should be instantiated by derived class see HydrogenPorousTransportLayer.m and OxygenPorousTransportLayer.m
            compInd.ncomp     = 5;
            compInd.liquid    = [compInd.H2Oliquid; compInd.OH; compInd.K];
            compInd.gas       = [compInd.H2Ogas; compInd.activeGas];
            
            phaseInd.liquid = 1;
            phaseInd.gas    = 2;
            phaseInd.solid  = 3;
            phaseInd.mobile = [1; 2];
            phaseInd.nphase = 3;

            mobPhaseInd.liquid    = 1;
            mobPhaseInd.gas       = 2;
            mobPhaseInd.nmobphase = 2;
            mobPhaseInd.phaseMap  = [1; 2];
            
            % compInd.phaseMap(compInd.H2Oliquid)  = phaseInd.liquid;
            compInd.phaseMap  = [1; 2; 1; 1; 2]; % first component (H2Oliquid) is in phase indexed by 1 (liquid phase), and so on
            
            liquidInd.H2Oliquid = 1;
            liquidInd.OH = 2;
            liquidInd.K  = 3;
            liquidInd.ncomp  = 3;
            liqudInd.compMap = [1; 3; 4];
            
            gasInd.H2Ogas    = 1;
            gasInd.activeGas = 2;
            gasInd.ncomp     = 2;
            gasInd.compMap   = [2; 5];

            model.compInd     = compInd;
            model.phaseInd    = phaseInd;
            model.mobPhaseInd = mobPhaseInd;
            model.liquidInd   = liquidInd;
            model.gasInd      = gasInd;

            paramobj.Boundary.compInd     = compInd;
            paramobj.Boundary.phaseInd    = phaseInd;
            paramobj.Boundary.mobPhaseInd = mobPhaseInd;
            paramobj.Boundary.liquidInd   = liquidInd;
            paramobj.Boundary.gasInd      = gasInd;
            
            model.Boundary = PorousTransportLayerBoundary(paramobj.Boundary);
        end

        function model = registerVarAndPropfuncNames(model)
            
        %% Declaration of the Dynamical Variables and Function of the model
        % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@ElectronicComponent(model);

            phaseInd    = model.phaseInd;
            liquidInd   = model.liquidInd;
            gasInd      = model.gasInd;
            compInd     = model.compInd;
            mobphaseInd = model.mobPhaseInd;
            
            ncomp   = compInd.ncomp;
            ngas    = gasInd.ncomp;
            nph     = phaseInd.nphase;
            nliquid = liquidInd.ncomp;
            nmobph  = mobphaseInd.nmobphase;

            varnames = {};
            % Total concentration of OH- (mol per total volume, that is not only liquid volume) [mol m^-3]
            varnames{end + 1} = 'OHceps';
            % Total density of gas H2O  (mass per total volume, that is not only gas volume) in [kg m^-3]
            varnames{end + 1} = 'H2Ogasrhoeps';
            % Liquid volume fraction, without unit [-]
            varnames{end + 1} = 'liqeps';
            % total liquid density (mass of liquid per total volume) in [kg m^-3]
            varnames{end + 1} = 'liqrhoeps';
            % Phase pressures in [Pa]
            phasePressures = VarName({}, 'phasePressures', nph);
            varnames{end + 1} = phasePressures;
            % Phase volume fractions, without unit [-]
            volumeFractions = VarName({}, 'volumeFractions', nph);
            varnames{end + 1} = volumeFractions;
            % Partial pressure for each component of gas in [Pa]
            varnames{end + 1}  = VarName({}, 'compGasPressures', ngas);
            % Masses for each component of gas (per total volume, as used in mass of conservation law) in [kg m^-3]
            compGasMasses = VarName({}, 'compGasMasses', ngas);
            varnames{end + 1} = compGasMasses;
            % Liquid density (Mass of liquid per volume of liquid) in [kg m^-3]
            varnames{end + 1} = 'liqrho';
            % Concentrations in the liquid in [mol m^-3]
            concentrations = VarName({}, 'concentrations', nliquid);
            varnames{end + 1} = concentrations;
            % Concentrations in the OH molality
            varnames{end + 1} = 'OHmolality';
            % H2O activity (we could assemble all activities there but it seems that only H2O activity is needed)
            varnames{end + 1} = 'H2Oa';
            % partialMolarVolumes = VarName({}, 'partialMolarVolumes', nliquid);
            % varnames{end + 1} = partialMolarVolumes;

            %% Flow variables

            viscosities = VarName({}, 'viscosities', nmobph);
            varnames{end + 1} = viscosities;            
            
            %% Phase velocities

            % phase velocity in [m s^-1] integrated over each cell face of the grid. Hence, unit is [m^3 s^-1]
            phaseFluxes = VarName({}, 'phaseFluxes', nph);
            varnames{end + 1} = phaseFluxes;
            
            %% Fluxes

            % Mass fluxes for the gass components in [kg s^-1] (integrated for each cell face in the grid)
            compGasFluxes = VarName({}, 'compGasFluxes', ngas);
            varnames{end + 1} = compGasFluxes;
            % Convective flux for OH in [mol s^-1] (unit is such because the flux is integrated for each cell face in the grid)
            varnames{end + 1} = 'convOHFlux';
            % Diffusion flux for OH in [mol s^-1] (unit is such integrated for each cell face in the grid)
            varnames{end + 1} = 'diffOHFlux';
            % Migration flux for OH in [mol s^-1] (unit is such integrated for each cell face in the grid)
            varnames{end + 1} = 'migOHFlux';
            % Mass flux for total of liquid components in [kg s^-1] (unit is such integrated for each cell face in the grid)
            varnames{end + 1} = 'liquidMassFlux';
            
            % Vapor pressure in [Pa]
            varnames{end + 1} = 'vaporPressure';
            
            
            %% Coupling variables
            
            % Mass sources for the gas Components in [kg s^-1] (source term for each grid cell)
            varnames{end + 1}  = VarName({}, 'compGasSources', ngas);
            % Mass sources at the boundaries for the gas components in [kg s^-1] (source term for each grid cell)
            varnames{end + 1}  = VarName({}, 'compGasBcSources', ngas);
            % Accumulation term for the gass components in [kg s^-1]
            varnames{end + 1}  = VarName({}, 'compGasAccums', ngas);
            % Source of OH in [mol s^-1]  (source term for each grid cell)
            varnames{end + 1} = 'OHSource';
            % Source of OH in [mol s^-1] at boundary  (yet one value for each grid cell)
            varnames{end + 1} = 'OHbcSource';
            % Mass Source of liquid in [kg s^-1] (source term for each grid cell)
            varnames{end + 1} = 'liquidSource';
            % Mass Source of liquid in [kg s^-1] at boundary (yet one value per for each grid cell))
            varnames{end + 1} = 'liquidBcSource';
            % Liquid-Vapor exchange rate for H2O (H2Ogas <-> H2Oliquid) in [mol m^-3 s^-1)]
            varnames{end + 1} = 'H2OvaporLiquidExchangeRate';
            % Source of H2Oliquid in [mol s^-1] (source term for each grid cell)
            varnames{end + 1} = 'H2OliquidSource';
            % Accumulation terms for OH in [mol s^-1] (accumulation term for each grid cell)
            varnames{end + 1} = 'OHaccum';
            % Accumulation term for the overall liquid components in [kg s^-1] (accumulation term for each grid cell)
            varnames{end + 1} = 'liquidAccumTerm';
            
            %% Residual variables
            
            % Residual for the mass conservation equations of the components in the gas
            varnames{end + 1}  = VarName({}, 'compGasMassCons', ngas);
            % Residual for the mass conservation equation of the aggregated components in the liquid phase
            varnames{end + 1} = 'liquidMassCons';
            % Residual for the conservation equation for OH in the liquid phase
            varnames{end + 1} = 'OHMassCons';
            % Residual for the equation of state of the liquid
            varnames{end + 1} = 'liquidStateEquation';

            model = model.registerVarNames(varnames);

            
            fn = @() PorousTransportLayer.updateVolumeFractions;
            inputnames = {'liqeps'};
            model = model.registerPropFunction({volumeFractions, fn, inputnames});            

            % assemble liquid pressure using capillary pressure function
            fn = @() PorousTransportLayer.updateLiquidPressure;
            inputnames = {VarName({}, 'phasePressures', nph, phaseInd.gas), volumeFractions};
            model = model.registerPropFunction({VarName({}, 'phasePressures', nph, phaseInd.liquid), fn, inputnames});
            
            %% assemble masses 
            
            % assemble mass of OH
            fn = @() PorousTransportLayer.updateOHconcentration;
            inputnames = {'OHceps', ...
                          VarName({}, 'volumeFractions', nph, phaseInd.liquid)};
            var = VarName({}, 'concentrations', nliquid, liquidInd.OH);
            model = model.registerPropFunction({var, fn, inputnames});
            
            % assemble mass of H2O in gas phase
            fn = @() PorousTransportLayer.updateMassH2Ogas;
            inputnames = {'H2Ogasrhoeps'};
            var = VarName({}, 'compGasMasses', ngas, gasInd.H2Ogas);
            model = model.registerPropFunction({var, fn, inputnames});            
            
            % update liquid density
            fn = @() PorousTransportLayer.updateLiquidDensity;
            inputnames = {'liqrhoeps', ...
                          VarName({}, 'volumeFractions', nph, phaseInd.liquid)};
            model = model.registerPropFunction({'liqrho', fn, inputnames});            

            % update conductivity
            fn = @() PorousTransportLayer.updateConductivity;
            inputnames = {'T', ...
                          VarName({}, 'concentrations', nliquid, liquidInd.OH)};
            model = model.registerPropFunction({'conductivity', fn, inputnames});            
            
            % assemble concentrations
            fn = @() PorousTransportLayer.updateConcentrations;
            inputnames = {'liqrho', ...
                          VarName({}, 'concentrations', nliquid, liquidInd.OH)};
            ind = setdiff([1 : nliquid]', liquidInd.OH);
            model = model.registerPropFunction({VarName({}, 'concentrations', nliquid, ind), fn, inputnames});
            
            fn = @() PorousTransportLayer.updateWaterActivity;
            inputnames = {'T', 'OHmolality'};
            model = model.registerPropFunction({'H2Oa', fn, inputnames});

            % compute OH molality
            fn = @() PorousTransportLayer.updateMolality;
            inputnames = {VarName({}, 'concentrations', nph, model.liquidInd.OH), 'liqrho'};
            model = model.registerPropFunction({'OHmolality', fn, inputnames});            
            
            % compute vapor pressure
            fn = @() PorousTransportLayer.updateVaporPressure;
            inputnames = {'T', 'OHmolality'};
            model = model.registerPropFunction({'vaporPressure', fn, inputnames});
            
            % update evaporation term
            fn = @() PorousTransportLayer.updateEvaporationTerm;
            inputnames = {'T', ...
                          'vaporPressure', ...
                          VarName({}, 'compGasPressures', ngas, gasInd.H2Ogas), ...
                          volumeFractions};
            model = model.registerPropFunction({'H2OvaporLiquidExchangeRate', fn, inputnames});

            
            % Assemble phase velocities
            fn = @() PorousTransportLayer.updatePhaseVelocities;
            for imobph = 1 : nmobph
                iphase = phaseInd.mobile(imobph);
                inputnames = {VarName({}, 'phasePressures', nph, iphase), ...
                              VarName({}, 'volumeFractions', nph, iphase), ...
                              VarName({}, 'viscosities', nmobph, imobph)};
                model = model.registerPropFunction({VarName({}, 'phaseFluxes', nph, iphase), fn, inputnames});
            end
            
            % assemble OH convection flux
            fn = @() PorousTransportLayer.updateOHConvectionFlux;
            inputnames = {VarName({}, 'concentrations', nliquid, liquidInd.OH), ...
                          VarName({}, 'phaseFluxes', nph, phaseInd.liquid)};
            model = model.registerPropFunction({'convOHFlux', fn, inputnames}); 
            
            % assemble OH diffusion flux
            fn = @() PorousTransportLayer.updateOHDiffusionFlux;
            inputnames = {VarName({}, 'concentrations', nliquid, liquidInd.OH), ...
                          VarName({}, 'volumeFractions', nph, phaseInd.liquid)};
            model = model.registerPropFunction({'diffOHFlux', fn, inputnames}); 
            
            % assemble OH migration flux
            fn = @() PorousTransportLayer.updateOHMigrationFlux;
            inputnames = {'j'};
            model = model.registerPropFunction({'migOHFlux', fn, inputnames}); 
            
            % assemble fluxes of gas components
            fn = @() PorousTransportLayer.updateGasFluxes;
            for igas = 1 : ngas
                inputnames = {VarName({}, 'compGasMasses', ngas, igas)         , ...
                              VarName({}, 'volumeFractions', nph, phaseInd.gas), ...
                              VarName({}, 'phaseFluxes', nph, phaseInd.gas)};
                model = model.registerPropFunction({VarName({}, 'compGasFluxes', ngas, igas), fn, inputnames});
            end
            
            % Assemble flux of the overall liquid component
            fn = @() PorousTransportLayer.updateLiquidMassFlux;
            inputnames = {VarName({}, 'phaseFluxes', nph, phaseInd.liquid), ...
                          'liqrho'};
            model = model.registerPropFunction({'liquidMassFlux', fn, inputnames});
                
            fn = @() PorousTransportLayer.updateESource;
            inputnames = {'OHSource'};
            model = model.registerPropFunction({'eSource', fn, inputnames});

            fn = @() PorousTransportLayer.updateLiquidSource;
            inputnames = {'OHSource', ...
                          'H2OliquidSource'};
            model = model.registerPropFunction({'liquidSource', fn, inputnames});

            fn = @() PorousTransportLayer.updateCurrent;
            inputnames = {'T'           , ....
                          'phi'         , ...
                          'conductivity', ...
                          VarName({}, 'concentrations', nliquid, liquidInd.OH)};
            model = model.registerPropFunction({'j', fn, inputnames});
            
            fn = @() PorousTransportLayer.updateLiquidViscosity;
            inputnames = {'T', ....
                          VarName({}, 'concentrations', nliquid, liquidInd.OH)};
            model = model.registerPropFunction({VarName({}, 'viscosities', nph, phaseInd.liquid), fn, inputnames});
            
            fn = @() PorousTransportLayer.updateH2OgasSource;
            inputnames = {'H2OvaporLiquidExchangeRate'};
            model = model.registerPropFunction({VarName({}, 'compGasSources', ngas, gasInd.H2Ogas), fn, inputnames});

            % Assemble mass conservation equations for components in gas phase
            fn = @() PorousTransportLayer.updateGasMassCons;
            for igas = 1 : ngas
                inputnames = {VarName({}, 'compGasBcSources', ngas, igas), ...
                              VarName({}, 'compGasFluxes'   , ngas, igas), ...
                              VarName({}, 'compGasSources'  , ngas, igas), ...
                              VarName({}, 'compGasAccums'   , ngas, igas)};
                model = model.registerPropFunction({VarName({}, 'compGasMassCons', ngas, igas), fn, inputnames});
            end

            % Assemble mass conservation for the overall liquid component
            fn = @() PorousTransportLayer.updateLiquidMassCons;
            inputnames = {'liquidMassFlux'     , ...
                          'liquidAccumTerm', ...
                          'liquidBcSource', ...
                          'liquidSource'};
            model = model.registerPropFunction({'liquidMassCons', fn, inputnames});
            
            % Assemble mass conservation equation for OH
            fn = @() PorousTransportLayer.updateOHMassCons;
            inputnames = {'convOHFlux', ...
                          'diffOHFlux', ...
                          'migOHFlux' , ...
                          'OHbcSource', ...
                          'OHSource'  , ...
                          'OHaccum'};
            model = model.registerPropFunction({'OHMassCons', fn, inputnames});

            % Assemble partial Molar Volumes (not used in first implementation)
            % fn = @() PorousTransportLayer.updatePartialMolarVolumes;
            % inputnames = {'OHmolality', 'T', 'liqrho'};
            % ind = [model.liquidInd.OH; model.liquidInd.K];
            % model = model.registerPropFunction({VarName({}, 'partialMolarVolumes', nliquid, ind), fn, inputnames});

            % we use liquid incompressibility for the moment
            fn = @() PorousTransportLayer.updateLiquidAccum;
            inputnames = {};
            model = model.registerPropFunction({'liquidAccumTerm', fn, inputnames});
            
            % Assemble residual of equation of state for the liquid phase
            fn = @() PorousTransportLayer.setupLiquidStateEquation;
            inputnames = {concentrations};
            model = model.registerPropFunction({'liquidStateEquation', fn, inputnames});

            fn = @() PorousTransportLayer.updateAccumTerms;
            functionCallSetupFn = @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction);
            fn = {fn, functionCallSetupFn};
            for igas = 1 : ngas
                inputnames = {VarName({}, 'compGasMasses', ngas, igas)};
                model = model.registerPropFunction({VarName({}, 'compGasAccums', ngas, igas), fn, inputnames});
            end
            inputnames = {'OHceps'};
            model = model.registerPropFunction({'OHaccum', fn, inputnames});

            % update BC terms
            bd = 'Boundary';

            fn = @() PorousTransportLayer.updateBcTerms;
            inputnames = {VarName({}, 'concentrations', nliquid, liquidInd.OH)                , ...
                          VarName({bd}, 'concentrations', nliquid, liquidInd.OH)              , ...
                          VarName({bd}, 'gasDensities', ngas)                                 , ...
                          'liqrho'                                                            , ...
                          {bd, 'liqrho'}                                                      , ...
                          'phi'                                                               , ...
                          {bd, 'phi'}                                                         , ...
                          VarName({}, 'viscosities', nmobph)                                  , ...
                          VarName({}, 'phasePressures', nph, [phaseInd.gas, phaseInd.liquid]) , ...
                          VarName({}, 'volumeFractions', nph, [phaseInd.gas, phaseInd.liquid]), ...
                          VarName({bd}, 'phasePressures', nph, [phaseInd.gas, phaseInd.liquid])
                         };
            model = model.registerPropFunction({VarName({}, 'compGasBcSources', ngas), fn, inputnames});
            model = model.registerPropFunction({'jBcSource', fn, inputnames});
            model = model.registerPropFunction({'OHbcSource', fn, inputnames});
            model = model.registerPropFunction({'liquidBcSource', fn, inputnames});
            model = model.registerPropFunction({VarName({bd}, 'bcEquations', 2 + ngas), fn, inputnames});
            
            model = model.removeVarName(VarName({}, 'phasePressures', nph, phaseInd.solid));
            model = model.removeVarName(VarName({}, 'phaseFluxes', nph, phaseInd.solid));
            
        end

        function state = updateBcTerms(model, state)
            
            coupterm   = model.externalCouplingTerm;
            harmFaceBc = model.operators.harmFaceBC;
            mphInd     = model.mobPhaseInd;
            phInd      = model.phaseInd;
            liqInd     = model.liquidInd;
            perm       = model.permeability;
            
            bcfaces = coupterm.couplingfaces;

            % Note: The fluxes that are computed here are all oriented towards interior domain (will come as source term in conservation equations)
            
            % We compute the fluxes of the mobile phases at the boundary

            for imph = 1 : mphInd.nmobphase
                iph = mphInd.phaseMap(imph);
                p    = state.phasePressures{iph};
                pBc  = state.(bd).phasePressures{iph};
                visc = state.viscosities{imph};
                vf   = state.volumeFractions{iph};
                [tc, bccells] = harmFaceBC((vf.^1.5).*perm./visc, bcfaces);
                phaseBcFlux{ph} = tc.*(pBc - p(bccells));
            end

            % We compute liquidBcSource (we use upwinding)

            liqrho   = state.liqrho;
            liqrhoBc = state.(bd).liqrho;
            
            liquibBcFlux = phaseBcFlux{phInd.liquid};
            upwindLiqrho = liqrho(bccells);
            upwindLiqrho(liquidBcFlux > 0) = liqrhoBc(liquidBcFlux > 0);

            state.liquidBcSource = upwindLiqrho*liquidBcFlux;

            % We compute electrical flux (should be consistent with what is implemented in method updateCurrent)

            phi   = state.phi;
            phiBc = state.(bd).phi;
            kappa = state.conductivity;
            tc = harmFaceBC(kappa, bcfaces);

            jBcFlux = tc.*(phiBc - phi(bccells));

            cOH = state.concentrations{liqInd.OH};
            cOHbc = state.cOH;
            
            dmudc = R.*T./cOH;
            R     = model.constants.R;
            F     = model.constants.F;
            t     = model.sp.OH.t;
            z     = model.sp.OH.z;             
            coef  = (conductivity./F).*(t/z).*dmudc;

            tc = harmFaceBC(coef, bcfaces);

            jBcFlux = jBcFlux + tc.*(cOHbc - cOH);
            
            % We compute the OH boundary fluxes
            
            % We compute OH boundary diffusion flux (should be consistent with what is implemented in method
            % updateOHDiffusionFlux)
            vf = state.volumeFractions{phInd.liquid};
            tc = harmFaceBC((vf.^1.5).*D, bcfaces);
            diffOHbcFlux = tc.*(cOHbc - cOH(bccells));
            
            % We compute OH boundary convection flux (should be consistent with what is implemented in method
            % updateOHConvectionFlux)
            upwindcOH = cOH;
            upwindcOHbc(liquidBcFlux > 0) = cOHbc(liquidBcFlux > 0);
            convOHbcFlux = upwindcOH.*liquiBcFlux;
            
            % We compute OH boundary migration flux (should be consistent with what is implemented in method
            % updateOHMigrationFlux)
            migOHbcFlux = t./(z.*F).*jBcFlux;

            OHbcSource = diffOHbcFlux + convOHbcFlux + migOHbcFlux;
            
            % We compute the boundary fluxes for the gases (should be consistent with what is implemented in method
            % updateGasFluxes)
            gasBcFlux = phaseBcFlux{phInd.gas};
            vf = state.volumeFractions{phInd.gas};
            for igas = 1 : numel(gasInd.ngas)
                rho = state.compGasMasses{igas}(bccells)./vf(bccells);
                rhoBc =  state.gasDensities{igas};
                upwindrhoBc = rho;
                upwindrhoBc(gasBcFlux > 0) = rhoBc(gasBcFlux > 0);
                compGasBcSources{igas} = rhoBc.*gasBcFlux;
            end


            % we update control
            for igas = 1 : ngas
                bceqs{igas}                = state.(bd).gasDensities{igas} - model.(bd).controlValues.gasDensities{igas};
                bceqs{igas}(gasBcFlux < 0) = state.(bd).gasDensities{igas} - state.gasDensities{igas}(bccells);
            end

            bceqs{igas + 1}                   = state.(bd).concentrations{liquidInd.OH} - model.(bd).controlValues.cOH;
            bceqs{igas + 1}(liquidBcFlux < 0) = state.(bd).concentrations{liquidInd.OH} - state.concentrations{liquidInd}(bccells);

            bceqs{igas + 2}                   = state.(bd).liqrho - model.(bd).liqrho;
            bceqs{igas + 2}(liquidBcFlux < 0) = state.(bd).liqrho - state.liqrho(bccells);
            

            % Each variable is first initialized (in a way we get AD )
            for igas = 1 : gasInd.ngas
                state.compGasBcSources{igas} = 0*phi;
                state.compGasBcSources{igas}(bccells) = compGasBcSources{igas};
            end
            state.OHbcSource = 0*phi;
            state.OHbcSource(bccells) = OHbcSource;
            state.jBcSource = 0*phi;
            state.jBcSource(bccells) = jBcSource;
            state.liquidBcSource = 0*phi;
            state.liquidBcSource(bccells) = liquidBcSource;
        end

        function state = updateVolumeFractions(model, state)

            liqeps = state.liqeps;

            state.volumeFractions{model.phaseInd.liquid} = liqeps;
            state.volumeFractions{model.phaseInd.solid}  = model.solidVolumeFraction;
            state.volumeFractions{model.phaseInd.gas}    = 1 - (liqeps + model.solidVolumeFraction);
            
        end

        function state = updateLiquidAccum(model, state)

        % We assume incompressibility for the moment
            state.liquidAccumTerm = 0;
            
        end

        function state = updateAccumTerms(model, state, state0, dt)

            vols   = model.G.cells.volumes;
            gasInd = model.gasInd;
            
            for igas = 1 : gasInd.ngas

                state.compGasAccums{igas} = vols.*(state.compGasMasses{igas} - state0.compGasMasses{igas})/dt;
                
            end
            
        end
        
        
        function state = updateLiquidViscosity(model, state)
        %   Data from Guo et al. Ref [2] Fitting parameters from Le
        %   Bideau, et al., Ref [3]. Valid from 20 to 60 degC and 2 to
        %   40 wt%. Units are kg m^-1 s^-1 (TBC), also known as Pa s.

            T   = state.T;
            cOH = state.concentrations{model.liquidInd.OH};
            
            coefs = [4.3e-1   ;
                     -2.51e-2 ;
                     10^-4    ;
                     1.3e-1 ];

            T   = T - 273.15; % convert to celcius
            cOH = 1e-3*cOH;   % convert to mol/litre
            
            mu = exp(coefs(1)       + ...
                     coefs(2).*T    + ...
                     coefs(3).*T.^2 + ...
                     coefs(4).*cOH);

            state.viscosities{model.mobPhaseInd.liquid} = mu;
            
        end
        
        function state = updateConductivity(model, state)
        %   Conductivity calculated according to the empirical model
        %   published by Gilliam, et al. Valid from 0 -
        %   12 M and 0 - 100 degC. Units are S m^-1.
            
            T   = state.T;
            cOH = state.concentrations{model.liquidInd.OH};
            
            coefs = [-2.041   ;
                     -0.0028  ;
                     0.005332 ;
                     207.2    ;
                     0.001043 ;
                     -0.0000003];

            cOH = 1e-3*cOH; % we convert to mol/litre s
            state.conductivity = (coefs(1).*cOH    + ...
                                  coefs(2).*cOH.^2 + ...
                                  coefs(3).*cOH.*T + ...
                                  coefs(4).*cOH./T + ...
                                  coefs(5).*cOH.^3 + ...
                                  coefs(6).*cOH.^2.*T.^2) .* 100;
            

        end
        
        function state = updateCurrent(model, state)

            R = model.constants.R;
            F = model.constants.F;
            t = model.sp.OH.t;
            z = model.sp.OH.z;             
            
            cOH   = state.concentrations{model.liquidInd.OH};
            kappa = state.conductivity;
            T     = state.T;
            phi   = state.phi;

            j = assembleFlux(phi, kappa);

            dmudc = R.*T./cOH;
            coef = (conductivity./F).*(t/z).*dmudc;

            j = j + assembleFlux(cOH, coef);

            state.j = j;
            
        end
        
        
        function state = updateLiquidDensity(model, state)

            vf = state.volumeFractions{model.phaseInd.liquid};

            state.liqrho = state.liqrhoeps./vf;
            
        end

        
        function state = updateLiquidPressure(model, state)
        % assemble liquid pressure using capillary pressure function

            K        = model.Permeability;
            levcoefs = model.leverettCoefficients;
            theta    = model.theta;
            
            pgas = state.phasePressures{model.phaseInd.gas};

            vl = state.volumeFractions{model.phaseInd.liquid};
            vg = state.volumeFractions{model.phaseInd.gas};
            vs = state.volumeFractions{model.phaseInd.solid};
            
            % Liquid saturation
            s = vl./(vl + vg);
            
            pc = 0.0694 .* cosd(theta) ./ sqrt(K./vs) .* leverett(levcoefs, s);
                        
            pliq = pgas - pc;

            state.phasePressures{model.phaseInd.liquid} = pliq;
        end

        
        function state = updateOHconcentration(model, state)

            vf = state.volumeFractions{model.phaseInd.liquid};
            OHceps = state.OHceps;
            
            state.concentrations{model.liquidInd.OH} = OHceps./vf;
        end
        
        
        function state = updateMassH2Ogas(model, state)
        % simple copy
            state.compGasMasses{compind.H2Ogas} = state.H2Ogasrhoeps;
            
        end

        function state = updateWaterActivity(model, state)

            T = state.T;
            m = state.OHmolality;

            % From Balej 1985 (ref 8), equation 28
            state.H2Oa = 10.^(-0.02255 .* m + 0.001434 .* m.^2 + (1.38.*m - 0.9254.*m.^2)./T);
            
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
            kLV   = model.sp.H2O.kLV;
            R     = model.constants.R;

            pH2Ovap = state.vaporPressure;
            T       = state.T;
            pH2Ogas = state.compGasPressures{model.gasInd.H2Ogas};
            vl      = state.volumeFractions{model.phaseInd.liquid};
            vg      = state.volumeFractions{model.phaseInd.gas};
            
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
        
        function state = updatePhaseVelocities(model, state)
        %% Assemble phase velocities

            K    = model.Permeability;
            pmap = model.mobPhaseInd.phaseMap;
            
            phaseinds = [model.phaseInd.liquid; model.phaseInd.gas];
            
            for ind = 1 : numel(phaseinds)
                phind = phaseinds{ind};
                p = state.phasePressures{phind};
                mu = state.viscosities{pmap(phind)};
                vf = state.volumeFractions{phins};
                v{ind} = assembleFlux(model, p, (vf.^1.5).*K./mu);
            end
            
            state.phaseFluxes = v;
            
        end
        
        function state = updateOHConvectionFlux(model, state)

            cOH = state.concentrations{model.liquidInd.OH};
            v   = state.phaseFluxes{model.phaseInd.liquid};

            state.convOHFlux = assembleUpwindFlux(model, v, cOH);
            
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
            v  = state.phaseFluxes{model.phaseInd.gas};

            for igas = 1 : model.gasInd.ncomp
                state.compGasFluxes{igas} =  assembleUpwindFlux(model, v, state.compGasMasses{igas}./vf);
            end
            
        end
        
        function state = updateLiquidMassFlux(model, state)
             
            rho = state.liqrho;
            v   = state.phaseFluxes{model.phaseInd.liquid};

            state.liquidMassFlux = assembleUpwindFlux(model, v, rho);
            
        end

        function state = updateMolality(model, state)
            
            MW = model.sp.OH.MW;
            
            rho = state.liqrho;
            cOH = state.concentrations{model.liquidInd.OH};
            
            state.OHmolality = cOH./(rho - cOH.*MW);
            
        end
        
        
        function state = updateESource(model, state)
        % Assemble charge source

            F = model.constants.F;
            z = model.sp.OH.z
            
            OHsrc = state.OHSource;
            
            state.eSource = F*z*OHsrc;
            
        end
        

        function state = updateGasMassCons(model, state)
        % Assemble mass conservation equations for components in gas phase
            
            for igas = 1 : model.gasInd.ncomp
                state.compGasMassCons{igas} = assembleConservationEquation(model, ...
                                                                           state.compGasFluxes{igas}   , ...
                                                                           state.compGasBcSources{igas}, ...
                                                                           state.compGasSources{igas}  , ...
                                                                           state.compGasAccums{igas});
            end
            
        end

        
        function state = updateGasMassCons0(model, state)
        % zero flux version
            for ind = 1 : model.gasInd.ncomp
                state.compGasFluxes{ind} = 0;
            end
            
            state = model.updateGasMassCons(state);
            
        end

        function state = updateLiquidSource(model, state)

            % Note that the units are
            % OHsource                   : [mol/s]
            % H2OliquidSource            : [mol/s]
            
            state.liquidSource = state.OHsource.*(model.sp.OH.MW + model.sp.K.MW) + state.H2OliquidSource.*model.sp.H2O.MW;

        end

        function state = updateH2OgasSource(model, state)

            vols = model.G.cells.volumes;
            
            state.compGasSources{model.gasInd.H2O} = -model.sp.H2O.MW*vols.*state.H2OvaporLiquidExchangeRate;

        end
        
        function state = updateLiquidMassCons(model, state)
        % Assemble mass conservation for the overall liquid component
            
            state.liquiMassCons = assembleConservationEquation(model             , ...
                                                               state.liquidMassFlux  , ...
                                                               state.liquidSource, ...
                                                               state.liquidBcSource, ...
                                                               state.liquidAccumTerm);
            
        end

        function state = updateOHMassCons(model, state)
        % Assemble mass conservation equation for OH

            OHflux = state.convOHFlux + state.diffOHFlux + state.migOHFlux;
            state.OHMassCons = assembleConservationEquation(model           , ...
                                                            OHFlux          , ...
                                                            state.OHbcSource, ...
                                                            state.OHsource  , ...
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


        function state = updateConcentrations(model, state)
            
            sp = model.sp;

            cOH = state.concentrations{model.liquidInd.OH};
            liqrho = state.liqrho;

            cK   = cOH;
            cH2O = (liqrho - cOH.*sp.OH.MW - cK.*sp.K.MW)./sp.H2O.MW;
            % cH   = 1e3.*(10.^-sp.H2O.beta .* (1e-3.*cOH).^-1);
            
            lInd = model.liquidInd;
            state.concentrations{lInd.K}   = cK;
            state.concentrations{lInd.H2O} = cH2O;
            % TODO  : check if we need H+ concentration later, if yes, we should uncomment the line below and adjust the indexings.
            % state.concentrations{lInd.H}   = cH;
            
        end

        function state = setupLiquidStateEquation(model, state)
            
            liqStateEq = -1;
            for ind = 1 : model.liquidInd.ncomp
                liqStateEq = liqStateEq + state.concentrations{ind}.*model.V0s(ind);
            end
            
            state.liquidStateEquation = liqStateEq;
            
        end
        
    end


end

