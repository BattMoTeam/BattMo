classdef PorousTransportLayer < ElectronicComponent

    properties


        Boundary

        compInd     % mapping structure for component indices
        phaseInd    % mapping structure for phase indices
        mobPhaseInd % mapping structure for mobile phase indices
        liquidInd   % mapping structure for component indices
        gasInd      % mapping structure for component indices

        solidVolumeFraction  % Solid volume fraction (should correspond to total ionomer + inactive part)
        leverettCoefficients % coefficients for the Leverett function definition, see method updateLiquidPressure and function leverett.m
        theta                % water contact angle
        permeability         % Permeability [Darcy]
        tortuosity           % Tortuosity

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

        externalCouplingTerm

        % Valuer for partial molar volumes of the liquid components
        % In the current implementation, the values are computed and set using the initial data.
        Vs

        % helper structures (those are not given as input but initialized with the model)
        MW    % molecular weight of KOH solution
        gasMW % Molecular weight of the active gas (H2 og O2). It will be initialized by the child class (HydrogenPorousTransportLayer or OxygenPorousTransportLayer).

    end

    methods

        function model = PorousTransportLayer(inputparams)

            % inputparams is instance of ElectronicComponentInputParams
            inputparams.use_thermal = false;
            model = model@ElectronicComponent(inputparams);

            %% Setup the model using the input parameters
            fdnames = {'G'                   , ...
                       'solidVolumeFraction' , ...
                       'leverettCoefficients', ...
                       'theta'               , ...
                       'tortuosity'          , ...
                       'permeability'        , ...
                       'sp'                  , ...
                       'MW'                  , ...
                       'externalCouplingTerm'};
            model = dispatchParams(model, inputparams, fdnames);


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

            % compInd.phaseMap(compInd.H2O)  = phaseInd.liquid;
            compInd.phaseMap  = [1; 2; 1; 1; 2]; % first component (H2O) is in phase indexed by 1 (liquid phase), and so on

            liquidInd.H2O     = 1;
            liquidInd.OH      = 2;
            liquidInd.K       = 3;
            liquidInd.nliquid = 3;
            liqudInd.compMap  = [1; 3; 4];

            gasInd.H2O       = 1;
            gasInd.activeGas = 2;
            gasInd.ngas      = 2;
            gasInd.compMap   = [2; 5];

            model.compInd     = compInd;
            model.phaseInd    = phaseInd;
            model.mobPhaseInd = mobPhaseInd;
            model.liquidInd   = liquidInd;
            model.gasInd      = gasInd;

            inputparams.Boundary.compInd     = compInd;
            inputparams.Boundary.phaseInd    = phaseInd;
            inputparams.Boundary.mobPhaseInd = mobPhaseInd;
            inputparams.Boundary.liquidInd   = liquidInd;
            inputparams.Boundary.gasInd      = gasInd;

            model.Boundary = PorousTransportLayerBoundary(inputparams.Boundary);

            % initialize helping structures
            sp = model.sp;
            model.MW = sp.OH.MW + sp.K.MW;

        end

        function model = registerVarAndPropfuncNames(model)

        %% Declaration of the Dynamical Variables and Function of the model
        % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@ElectronicComponent(model);

            phaseInd    = model.phaseInd;
            liquidInd   = model.liquidInd;
            gasInd      = model.gasInd;
            compInd     = model.compInd;
            mobPhaseInd = model.mobPhaseInd;

            ncomp   = compInd.ncomp;
            ngas    = gasInd.ngas;
            nph     = phaseInd.nphase;
            nliquid = liquidInd.nliquid;
            nmobph  = mobPhaseInd.nmobphase;

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
            % potential in the porous transport layer (for the moment a constant, we assume infinite conductivity in PTL)
            varnames{end + 1} = 'E';
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
            varnames{end + 1} = 'OHsource';
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
            var = VarName({}, 'compGasMasses', ngas, gasInd.H2O);
            model = model.registerPropFunction({var, fn, inputnames});

            % update liquid density
            fn = @() PorousTransportLayer.updateLiquidDensity;
            inputnames = {'liqrhoeps', ...
                          VarName({}, 'volumeFractions', nph, phaseInd.liquid)};
            model = model.registerPropFunction({'liqrho', fn, inputnames});

            % update conductivity
            fn = @() PorousTransportLayer.updateConductivity;
            inputnames = {'T', ...
                          VarName({}, 'concentrations', nliquid, liquidInd.OH), ...
                          VarName({}, 'volumeFractions', nph, phaseInd.liquid)};
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
                          VarName({}, 'compGasPressures', ngas, gasInd.H2O), ...
                          volumeFractions};
            % evaporation rate H2O(l) <->> H2O(g)
            % Here, the sign indicated by the repeated arrow sign corresponds to positive sign of rate
            model = model.registerPropFunction({'H2OvaporLiquidExchangeRate', fn, inputnames});


            % Assemble phase velocities
            fn = @() PorousTransportLayer.updatePhaseVelocities;
            phaseinds = [model.phaseInd.liquid; model.phaseInd.gas];
            inputnames = {VarName({}, 'phasePressures', nph, phaseinds), ...
                          VarName({}, 'volumeFractions', nph, phaseinds), ...
                          VarName({}, 'viscosities', nmobph)};
            model = model.registerPropFunction({VarName({}, 'phaseFluxes', nph, phaseinds), fn, inputnames});

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
            inputnames = {VarName({}, 'compGasMasses', ngas)         , ...
                          VarName({}, 'volumeFractions', nph, phaseInd.gas), ...
                          VarName({}, 'phaseFluxes', nph, phaseInd.gas)};
            model = model.registerPropFunction({VarName({}, 'compGasFluxes', ngas), fn, inputnames});

            % Assemble flux of the overall liquid component
            fn = @() PorousTransportLayer.updateLiquidMassFlux;
            inputnames = {VarName({}, 'phaseFluxes', nph, phaseInd.liquid), ...
                          'liqrho'};
            model = model.registerPropFunction({'liquidMassFlux', fn, inputnames});

            fn = @() PorousTransportLayer.updateESource;
            inputnames = {'OHsource'};
            model = model.registerPropFunction({'eSource', fn, inputnames});

            fn = @() PorousTransportLayer.updateLiquidSource;
            inputnames = {'OHsource', ...
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
            model = model.registerPropFunction({VarName({}, 'compGasSources', ngas, gasInd.H2O), fn, inputnames});

            % Assemble mass conservation equations for components in gas phase
            fn = @() PorousTransportLayer.updateGasMassCons;
            inputnames = {VarName({}, 'compGasBcSources', ngas), ...
                          VarName({}, 'compGasFluxes'   , ngas), ...
                          VarName({}, 'compGasSources'  , ngas), ...
                          VarName({}, 'compGasAccums'   , ngas)};
            model = model.registerPropFunction({VarName({}, 'compGasMassCons', ngas), fn, inputnames});

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
                          'OHsource'  , ...
                          'OHaccum'};
            model = model.registerPropFunction({'OHMassCons', fn, inputnames});

            % Assemble partial Molar Volumes (not used in first implementation)
            % fn = @() PorousTransportLayer.updatePartialMolarVolumes;
            % inputnames = {'OHmolality', 'T', 'liqrho'};
            % ind = [model.liquidInd.OH; model.liquidInd.K];
            % model = model.registerPropFunction({VarName({}, 'partialMolarVolumes', nliquid, ind), fn, inputnames});

            fn = @() PorousTransportLayer.updateLiquidAccum;
            functionCallSetupFn = @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction);
            fn = {fn, functionCallSetupFn};
            inputnames = {'liqrhoeps'};
            model = model.registerPropFunction({'liquidAccumTerm', fn, inputnames});

            % Assemble residual of equation of state for the liquid phase
            fn = @() PorousTransportLayer.setupLiquidStateEquation;
            inputnames = {concentrations};
            model = model.registerPropFunction({'liquidStateEquation', fn, inputnames});

            fn = @() PorousTransportLayer.updateAccumTerms;
            functionCallSetupFn = @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction);
            fn = {fn, functionCallSetupFn};
            inputnames = {VarName({}, 'compGasMasses', ngas), 'OHceps'};
            model = model.registerPropFunction({VarName({}, 'compGasAccums', ngas), fn, inputnames});
            model = model.registerPropFunction({'OHaccum', fn, inputnames});

            % update BC terms
            bd = 'Boundary';

            fn = @() PorousTransportLayer.updateBcTerms;
            inputnames = {VarName({}, 'concentrations', nliquid, liquidInd.OH)                , ...
                          {bd, 'cOH'}                                                               , ...
                          VarName({bd}, 'gasDensities', ngas)                                 , ...
                          'liqrho'                                                            , ...
                          {bd, 'liqrho'}                                                      , ...
                          VarName({}, 'viscosities', nmobph)                                  , ...
                          VarName({}, 'phasePressures', nph, [phaseInd.gas, phaseInd.liquid]) , ...
                          VarName({}, 'volumeFractions', nph, [phaseInd.gas, phaseInd.liquid]), ...
                          VarName({bd}, 'phasePressures', nph, [phaseInd.gas, phaseInd.liquid])
                         };
            model = model.registerPropFunction({VarName({}, 'compGasBcSources', ngas), fn, inputnames});
            model = model.registerPropFunction({'OHbcSource', fn, inputnames});
            model = model.registerPropFunction({'liquidBcSource', fn, inputnames});
            model = model.registerPropFunction({VarName({bd}, 'bcEquations', 2 + ngas), fn, inputnames});

            model = model.removeVarName(VarName({}, 'phasePressures', nph, phaseInd.solid));
            model = model.removeVarName(VarName({}, 'phaseFluxes', nph, phaseInd.solid));

        end

        function state = updateBcTerms(model, state)

            bd = 'Boundary';

            coupterm   = model.externalCouplingTerm;
            mphInd     = model.mobPhaseInd;
            phInd      = model.phaseInd;
            liquidInd  = model.liquidInd;
            gasInd     = model.gasInd;
            perm       = model.permeability;
            tau        = model.tortuosity;

            bcfaces = coupterm.couplingfaces;

            % Note: The fluxes that are computed here are all oriented towards interior domain (will come as source term in conservation equations)

            % We compute the fluxes of the mobile phases at the boundary

            for imph = 1 : mphInd.nmobphase
                iph = mphInd.phaseMap(imph);
                p    = state.phasePressures{iph};
                pBc  = state.(bd).phasePressures{iph};
                visc = state.viscosities{imph};
                vf   = state.volumeFractions{iph};
                [tc, bccells] = model.G.getBcHarmFace((vf.^tau).*perm./visc, bcfaces);
                phaseBcFlux{iph} = tc.*(pBc - p(bccells));
            end

            % We compute liquidBcSource (we use upwinding)

            liqrho   = state.liqrho;
            liqrhoBc = state.(bd).liqrho;

            liquidBcFlux = phaseBcFlux{phInd.liquid};

            % upwind version:
            % upwindFlag = liquidBcFlux > 0;
            % upwindLiqrho             = liqrho(bccells);
            % upwindLiqrho(upwindFlag) = liqrhoBc(upwindFlag);
            % liquidBcSource = upwindLiqrho*liquidBcFlux;

            % non-upwind version:
            liqrho = liqrho(bccells);
            liquidBcSource = liqrho*liquidBcFlux;

            %%%%%%%%%%%%%
            %
            % Expressions for electrical flux at boundary (should be consistent with what is implemented in method updateCurrent)
            % We do not use them because we assume zero electrical flux at boundary.
            % phi   = state.phi;
            % phiBc = state.(bd).phi;
            % kappa = state.conductivity;
            % tc = op.harmFaceBC(kappa, bcfaces);

            % jBcFlux = tc.*(phiBc - phi(bccells));

            % cOH = state.concentrations{liquidInd.OH};
            % cOHbc = state.cOH;

            % dmudc = R.*T./cOH;
            % R     = model.constants.R;
            % F     = model.constants.F;
            % t     = model.sp.OH.t;
            % z     = model.sp.OH.z;
            % coef  = (conductivity./F).*(t/z).*dmudc;

            % tc = op.harmFaceBC(coef, bcfaces);

            % jBcFlux = jBcFlux + tc.*(cOHbc - cOH);
            %
            %%%%%%%%%%%%%%

            % We compute the OH boundary fluxes

            % We compute OH boundary diffusion flux (should be consistent with what is implemented in method
            % updateOHDiffusionFlux)
            cOH = state.concentrations{liquidInd.OH};
            cOHbc = state.(bd).cOH;

            vf = state.volumeFractions{phInd.liquid};
            D  = model.sp.OH.D;
            tc = op.harmFaceBC((vf.^tau).*D, bcfaces);
            diffOHbcFlux = tc.*(cOHbc - cOH(bccells));

            % We compute OH boundary convection flux (should be consistent with what is implemented in method
            % updateOHConvectionFlux)

            % upwind version
            % upwindFlag = liquidBcFlux > 0;
            % upwindcOH             = cOH(bccells);
            % upwindcOH(upwindFlag) = cOHbc(upwindFlag);
            % convOHbcFlux = upwindcOH.*liquidBcFlux;
            % non-upwind version
            convOHbcFlux = cOH(bccells).*liquidBcFlux;

            % The OH boundary migration flux (should be consistent with what is implemented in method
            % updateOHMigrationFlux) is given by: migOHbcFlux = t./(z.*F).*jBcFlux
            % Since we assume no current at boundary it is equal to zero.
            migOHbcFlux = 0;

            OHbcSource = diffOHbcFlux + convOHbcFlux + migOHbcFlux;

            % We compute the boundary fluxes for the gases (should be consistent with what is implemented in method
            % updateGasFluxes)
            gasBcFlux = phaseBcFlux{phInd.gas};
            vf = state.volumeFractions{phInd.gas};
            for igas = 1 : gasInd.ngas

                rho   = state.compGasMasses{igas}./vf;
                rhoBc = state.(bd).gasDensities{igas};

                % upwind version
                % upwindFlag = gasBcFlux > 0;
                % upwindrho             = rho(bccells);
                % upwindrho(upwindFlag) = rhoBc(upwindFlag);
                % compGasBcSources{igas} = upwindrho.*gasBcFlux;

                % non-upwind version
                compGasBcSources{igas} = rho(bccells).*gasBcFlux;

                gasDensities{igas} = rho(bccells); % needed below

            end

            % we update control
            for igas = 1 : gasInd.ngas

                % upwind version
                % bcDirFlag = gasBcFlux > 0;
                % bcEquations{igas}            = state.(bd).gasDensities{igas} - gasDensities{igas};
                % bcEquations{igas}(bcDirFlag) = state.(bd).gasDensities{igas}(bcDirFlag) - model.(bd).controlValues.gasDensities{igas}(bcDirFlag);

                % non-upwind version
                bcEquations{igas} = state.(bd).gasDensities{igas} - gasDensities{igas};

            end

            bcEquations{igas + 1} = state.(bd).cOH - model.(bd).controlValues.cOH;


            % upwind version
            % bcDirFlag = liquidBcFlux > 0;
            % bcEquations{igas + 2}            = state.(bd).liqrho - state.liqrho(bccells);
            % bcEquations{igas + 2}(bcDirFlag) = state.(bd).liqrho(bcDirFlag) - model.(bd).controlValues.liqrho(bcDirFlag);

            % non-upwind version
            bcEquations{igas + 2} = state.(bd).liqrho - state.liqrho(bccells);

            % Each variable is first initialized (in a way we get AD )
            phi = state.phi;
            for igas = 1 : gasInd.ngas
                state.compGasBcSources{igas}          = 0*phi;
                state.compGasBcSources{igas}(bccells) = compGasBcSources{igas};
            end
            state.OHbcSource              = 0*phi;
            state.OHbcSource(bccells)     = OHbcSource;
            state.liquidBcSource          = 0*phi;
            state.liquidBcSource(bccells) = liquidBcSource;

            state.(bd).bcEquations = bcEquations;

        end

        function state = updateVolumeFractions(model, state)

            liqeps = state.liqeps;

            state.volumeFractions{model.phaseInd.liquid} = liqeps;
            state.volumeFractions{model.phaseInd.solid}  = model.solidVolumeFraction;
            state.volumeFractions{model.phaseInd.gas}    = 1 - (liqeps + model.solidVolumeFraction);

        end

        function state = updateLiquidAccum(model, state, state0, dt)

            vols   = model.G.getVolumes();

            state.liquidAccumTerm = vols.*(state.liqrhoeps - state0.liqrhoeps)/dt;

        end

        function state = updateAccumTerms(model, state, state0, dt)

            vols   = model.G.getVolumes();
            gasInd = model.gasInd;

            for igas = 1 : gasInd.ngas

                state.compGasAccums{igas} = vols.*(state.compGasMasses{igas} - state0.compGasMasses{igas})/dt;

            end

            state.OHaccum = vols.*(state.OHceps - state0.OHceps)/dt;

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

            mu = 1e-3*exp(coefs(1) + coefs(2)*T + coefs(3)*T.^2 + coefs(4)*cOH);

            state.viscosities{model.mobPhaseInd.liquid} = mu;

        end

        function state = updateConductivity(model, state)
        %   Conductivity calculated according to the empirical model
        %   published by Gilliam, et al. Valid from 0 -
        %   12 M and 0 - 100 degC. Units are S m^-1.

            tau = model.tortuosity;

            T   = state.T;
            cOH = state.concentrations{model.liquidInd.OH};
            vf  = state.volumeFractions{model.phaseInd.liquid};

            coefs = [-2.041   ;
                     -0.0028  ;
                     0.005332 ;
                     207.2    ;
                     0.001043 ;
                     -0.0000003];

            cOH = 1e-3*cOH; % we convert to mol/litre

            kappa = (coefs(1).*cOH    + ...
                     coefs(2).*cOH.^2 + ...
                     coefs(3).*cOH.*T + ...
                     coefs(4).*cOH./T + ...
                     coefs(5).*cOH.^3 + ...
                     coefs(6).*cOH.^2.*T.^2) .* 100;

            state.conductivity = kappa.*(vf.^tau);

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

            j = assembleFlux(model, phi, kappa);

            dmudc = R.*T./cOH;
            coef = (kappa./F).*(t/z).*dmudc;

            j = j + assembleFlux(model, cOH, coef);

            state.j = j;

        end


        function state = updateLiquidDensity(model, state)

            vf = state.volumeFractions{model.phaseInd.liquid};

            state.liqrho = state.liqrhoeps./vf;

        end

        function state = updateLiquidPressure(model, state)
        % assemble liquid pressure using capillary pressure function

            K        = model.permeability;
            levcoefs = model.leverettCoefficients;
            theta    = model.theta;

            pgas = state.phasePressures{model.phaseInd.gas};

            vfl = state.volumeFractions{model.phaseInd.liquid};
            vfg = state.volumeFractions{model.phaseInd.gas};
            vfs = state.volumeFractions{model.phaseInd.solid};

            % Liquid saturation
            sLiq = vfl./(vfl + vfg);

            pc = 0.0694 .* cosd(theta) ./ sqrt(K./(1 - vfs)) .* leverett(levcoefs, sLiq);

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
            state.compGasMasses{model.gasInd.H2O} = state.H2Ogasrhoeps;

        end

        function state = updateWaterActivity(model, state)

            T = state.T;
            m = state.OHmolality;

            % From Balej 1985 (ref 8), equation 28
            state.H2Oa = 10.^(-0.02255.*m + 0.001434.*m.^2 + (1.38.*m - 0.9254.*m.^2)./T);
            % state.H2Oa = 1 + 0*m;

        end

        function state = updateVaporPressure(model, state)

            m = state.OHmolality;
            T = state.T;

            a = -1.508e-2    .* m ...
                -1.6788e-3   .* m.^2 ...
                + 2.25887e-5 .* m.^3;

            b = 1 ...
                - 1.2062e-3 .* m ...
                + 5.6024e-4 .* m.^2 ...
                - 7.8228e-6 .* m.^3;

            pw0 = 10.^(35.4462 - (3343.93./T) - 10.9.*log10(T) + 4.1645e-3.*T);

            state.vaporPressure = 10.^(a + b.*log10(pw0)).*1e5;
        end

        function state = updateEvaporationTerm(model, state)

            %% assemble evaporation term
            psat = model.sp.H2O;
            MW   = model.sp.H2O.MW;
            kLV  = model.sp.H2O.kLV;
            R    = model.constants.R;

            pH2Ovap = state.vaporPressure;
            T       = state.T;
            pH2Ogas = state.compGasPressures{model.gasInd.H2O};
            vfl      = state.volumeFractions{model.phaseInd.liquid};
            vfg      = state.volumeFractions{model.phaseInd.gas};

            % (H2O)_liquid <->> (H2O)_vapor
            % Here, the direction of the reaction that is indicated by the repeated arrow symbol corresponds to a positive computed reaction rate

            sLiq = vfl./(vfl + vfg);

            evapSrc = 0*pH2Ogas; % initialize evapSrc (useful when AD)
            ind = (pH2Ogas > pH2Ovap);
            if any(ind)
                evapSrc(ind) = (pH2Ovap(ind) - pH2Ogas(ind)) .* kLV ./MW;
            end
            if any(~ind)
                evapSrc(~ind) = sLiq(~ind) .* (pH2Ovap(~ind) - pH2Ogas(~ind)) .* kLV ./ MW;
            end

            state.H2OvaporLiquidExchangeRate = evapSrc;

        end

        function state = updatePhaseVelocities(model, state)
        %% Assemble phase velocities

            K    = model.permeability;
            pmap = model.mobPhaseInd.phaseMap;
            tau  = model.tortuosity;

            phaseInds = [model.phaseInd.liquid; model.phaseInd.gas];

            for ind = 1 : numel(phaseInds)
                phind = phaseInds(ind);
                p = state.phasePressures{phind};
                mu = state.viscosities{pmap(phind)};
                vf = state.volumeFractions{phind};
                v{phind} = assembleFlux(model, p, (vf.^tau.*K./mu));
            end

            state.phaseFluxes = v;

        end

        function state = updateOHConvectionFlux(model, state)


            G = model.G;
            
            cOH = state.concentrations{model.liquidInd.OH};
            v   = state.phaseFluxes{model.phaseInd.liquid};

            % upwind version
            % state.convOHFlux = assembleUpwindFlux(model, v, cOH);

            % non-upwind version
            state.convOHFlux = (1./G.getTrans()).*G.getHarmFace(cOH).*v;

        end

        function state = updateOHDiffusionFlux(model, state)
        % assemble OH diffusion flux

            D   = model.sp.OH.D;
            tau = model.tortuosity;

            cOH = state.concentrations{model.liquidInd.OH};
            vf = state.volumeFractions{model.phaseInd.liquid};

            state.diffOHFlux = assembleFlux(model, cOH, vf.^tau.*D);

        end

        function state = updateOHMigrationFlux(model, state)
        % assemble OH migration flux
            F = model.constants.F;
            t = model.sp.OH.t;
            z = model.sp.OH.z;

            j = state.j;

            state.migOHFlux = t./(z.*F).*j;
        end

        function state = updateGasFluxes(model, state)
        % assemble convection fluxes

            G = model.G;
            
            vf = state.volumeFractions{model.phaseInd.gas};
            v  = state.phaseFluxes{model.phaseInd.gas};

            for igas = 1 : model.gasInd.ngas
                % upwind version
                % state.compGasFluxes{igas} =  assembleUpwindFlux(model, v, state.compGasMasses{igas}./vf);
                % non-upwind version :
                state.compGasFluxes{igas} = (1./G.getTrans()).*G.getHarmFace(state.compGasMasses{igas}./vf).*v;
            end

        end

        function state = updateLiquidMassFlux(model, state)

            rho = state.liqrho;
            v   = state.phaseFluxes{model.phaseInd.liquid};

            % upwind version
            % state.liquidMassFlux = assembleUpwindFlux(model, v, rho);
            % non upwind version:
            op = model.operators;
            state.liquidMassFlux = (1./op.T).*op.harmFace(rho).*v;

        end

        function state = updateMolality(model, state)

            MW = model.MW;

            rho = state.liqrho;
            cOH = state.concentrations{model.liquidInd.OH};

            state.OHmolality = cOH./(rho - cOH.*MW);

        end


        function state = updateESource(model, state)
        % Assemble charge source

            F = model.constants.F;
            z = model.sp.OH.z;

            OHsrc = state.OHsource;

            state.eSource = F*z*OHsrc;

        end


        function state = updateGasMassCons(model, state)
        % Assemble mass conservation equations for components in gas phase

            for igas = 1 : model.gasInd.ngas
                state.compGasMassCons{igas} = assembleConservationEquation(model, ...
                                                                           state.compGasFluxes{igas}   , ...
                                                                           state.compGasBcSources{igas}, ...
                                                                           state.compGasSources{igas}  , ...
                                                                           state.compGasAccums{igas});
            end

        end

        function state = updateLiquidSource(model, state)

            % Note that the units are
            % OHsource                   : [mol/s]
            % H2OliquidSource            : [mol/s]

            state.liquidSource = state.OHsource.*(model.sp.OH.MW + model.sp.K.MW) + state.H2OliquidSource.*model.sp.H2O.MW;

        end

        function state = updateH2OgasSource(model, state)

            vols = model.G.getVolumes();

            state.compGasSources{model.gasInd.H2O} = model.sp.H2O.MW*vols.*state.H2OvaporLiquidExchangeRate;

        end

        function state = updateLiquidMassCons(model, state)
        % Assemble mass conservation for the overall liquid component

            state.liquidMassCons = assembleConservationEquation(model             , ...
                                                                state.liquidMassFlux  , ...
                                                                state.liquidSource, ...
                                                                state.liquidBcSource, ...
                                                                state.liquidAccumTerm);

        end

        function state = updateOHMassCons(model, state)
        % Assemble mass conservation equation for OH

            OHflux = state.convOHFlux + state.diffOHFlux + state.migOHFlux;
            state.OHMassCons = assembleConservationEquation(model           , ...
                                                            OHflux          , ...
                                                            state.OHbcSource, ...
                                                            state.OHsource  , ...
                                                            state.OHaccum);
        end


        function state = updatePartialMolarVolumes(model, state)

            error('Not used in the current implementation');

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

            liquidInd = model.liquidInd;
            state.concentrations{liquidInd.K}   = cK;
            state.concentrations{liquidInd.H2O} = cH2O;
            % TODO  : check if we need H+ concentration later, if yes, we should uncomment the line below and adjust the indexings.
            % state.concentrations{lInd.H}   = cH;

        end

        function state = setupLiquidStateEquation(model, state)

            liqStateEq = -1;
            for ind = 1 : model.liquidInd.nliquid
                liqStateEq = liqStateEq + state.concentrations{ind}.*model.Vs(ind);
            end

            state.liquidStateEquation = liqStateEq;

        end

        function rho = density(model, c, T)
        % NOTE :
        % In the current implementatiom this method is only used for initialization

        % Calculates the density of aqueous KOH solution as a
        % function of concentration (c) in mol/m3 and temperature (T) in
        % K.
        %
        % Density data is compiled from Zatysev et al., Ref [1], and
        % fit using the MATLAB curve fitting toolbox. Units are kg/m3
        % or g/L. Valid from 0 - 50 wt% and 0 - 100 degC.


            coefs = [ 794.7015;
                      0.0456546;
                      1.6355;
                      -8.392e-7;
                      1.659e-5;
                      -0.003205;
                      1.73e-11;
                      -9.7647e-12;
                      -3.3927e-08 ];

            % Empirical model for mass density, rho = f(c,T)
            rho = coefs(1)                                         + ...
                  coefs(2).*c + coefs(3).*T                        + ...
                  coefs(4).*c.^2 + coefs(5).*c.*T + coefs(6).*T.^2 + ...
                  coefs(7).*c.^3                                   + ...
                  coefs(8).*c.^2.*T + coefs(9).*c.*T.^2;

        end

        function [Vs, cH2O] = partialMolarVolume(model, c, rho, T)
        % NOTE :
        % In the current implementatiom this method is only used for initialization
        %
        % rho : liquid density [kg m^-3]
        % c : concentration [mol m^-3]

            MW = 0.0561056;  % Molecular Weight of KOH [kg mol^-1]

            sp = model.sp;
            lind = model.liquidInd;

            % Calculate molal concentration
            m =  c./(rho - c*MW);

            % introduced alias just for clarity in the expressions below
            cOH = c;
            cK  = c;
            mOH = m;
            mK  = m;

            coefs = [778.6106;
                     42.8554;
                     1.7400;
                     -1.4428;
                     0.04585;
                     -0.0034;
                     0.0197;
                     4.1087e-4;
                     -1.1139e-4];

            % Empirical model for first derivative of mass density with
            % respect to KOH molality

            drhodm = coefs(2)                     + ...
                     2.*coefs(4).*m + coefs(5).*T + ...
                     3.*coefs(7).*m.^2            + ...
                     2.*coefs(8).*m.*T + coefs(9).*T.^2;

            pmv = -1./(rho.^2).*drhodm.*(1 + m.*MW) + (1./rho).*MW;

            Vs(lind.OH)  = pmv.*abs(sp.OH.V0)./(abs(sp.OH.V0) + abs(sp.K.V0));
            Vs(lind.K)   = pmv.*abs(sp.K.V0)./(abs(sp.OH.V0) + abs(sp.K.V0));
            Vs(lind.H2O) = sp.H2O.MW.* (1./rho +...
                                        mOH.*(sp.OH.MW./rho - Vs(lind.OH)) +...
                                        mK.*(sp.K.MW ./rho - Vs(lind.K)));

            cH2O = (1 - cOH.*Vs(lind.OH) - cK.*Vs(lind.K))./Vs(lind.H2O);

        end

    end

end






%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
