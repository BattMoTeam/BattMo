classdef Electrolyser < BaseModel
    
    properties
        
        con = PhysicalConstants();

        % Components
        IonomerMembrane
        OxygenEvolutionElectrode
        HydrogenEvolutionElectrode
        
        couplingTerms
        couplingNames

        primaryVarNames
        funcCallList

    end
    
    methods
        
        function model = Electrolyser(paramobj)

            model = model@BaseModel();

            fdnames = {'G' , ...
                       'couplingTerms'};
            model = dispatchParams(model, paramobj, fdnames);

            model.OxygenEvolutionElectrode   = EvolutionElectrode(paramobj.OxygenEvolutionElectrode);
            model.HydrogenEvolutionElectrode = EvolutionElectrode(paramobj.HydrogenEvolutionElectrode);
            model.IonomerMembrane            = IonomerMembrane(paramobj.IonomerMembrane);

            % setup couplingNames
            model.couplingNames = cellfun(@(x) x.name, model.couplingTerms, 'uniformoutput', false);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {'T', ...
                        VarName({}, 'controlEqs', 2)};

            model = model.registerVarNames(varnames);

            inm = 'IonomerMembrane';
            her = 'HydrogenEvolutionElectrode';
            oer = 'OxygenEvolutionElectrode';

            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';
            ptl = 'PorousTransportLayer';
            
            fn = @Electrolyser.dispatchTemperature;
            inputvarnames = {'T'};
            model = model.registerPropFunction({{her, 'T'}, fn, inputvarnames});            
            model = model.registerPropFunction({{oer, 'T'}, fn, inputvarnames});
            model = model.registerPropFunction({{inm, 'T'}, fn, inputvarnames});
            
            fn = @Electrolyser.updateC;

            eldes = {her, oer};
            layers = {ctl, exl};
            
            fn = @Electrolyser.updateIonomerSources;
            inputvarnames = {};
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                inputvarnames{end + 1} = {elde, exl, 'H2OexchangeRate'};
                inputvarnames{end + 1} = {elde, ctl, 'inmrOHsource'};
                inputvarnames{end + 1} = {elde, exl, 'OHexchangeRate'};
            end
            model = model.registerPropFunction({{inm, 'H2OSource'}, fn, inputvarnames});
            model = model.registerPropFunction({{inm, 'OHsource'}, fn, inputvarnames});


            fn = @Electrolyser.dispatchIonomerToReactionLayers;
            inputvarnames = {{inm, 'phi'}, ...
                             {inm, 'cOH'}, ...
                             {inm, 'H2Oa'}};
            for ielde = 1 : numel(eldes)
                for ilayer = 1 : numel(layers)
                    elde = eldes{ielde};
                    layer = layers{ilayer};
                    model = model.registerPropFunction({{elde, layer, 'cOHinmr'}, fn, inputvarnames});
                    model = model.registerPropFunction({{elde, layer, 'phiInmr'}, fn, inputvarnames});
                    model = model.registerPropFunction({{elde, layer, 'H2OaInmr'}, fn, inputvarnames});
                end
            end

            fn = @Electrolyser.setupControl;
            fn = {fn, @PropFunction.drivingForceFuncCallSetupFn};
            inputvarnames = {{oer, ctl, 'I'}, ...
                             {her, ctl, 'E'}
                            };
            model = model.registerPropFunction({VarName({}, 'controlEqs', 2), fn, inputvarnames});

            model = model.registerStaticVarNames({{inm, 'jBcSource'}     , ...
                                                  {oer, ptl, 'jBcSource'}, ...
                                                  {her, ptl, 'jBcSource'}, ...
                                                  'T'});

            model = model.removeVarNames({{her, ctl, 'I'}, ...
                                          {her, ctl, 'eSource'}});
        end


        function model = validateModel(model, varargin)

            model = validateModel@BaseModel(model, varargin{:});
            model = model.setupComputationalGraph();
            
            cgt = model.computationalGraph;
            model.primaryVarNames = cgt.getPrimaryVariables();
            model.funcCallList = cgt.setOrderedFunctionCallList();
            
        end
        
        function [model, state] = setupBcAndInitialState(model)
            
            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';
            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            inm = 'IonomerMembrane';
            
            con = model.con;

            pGas = 101325*Pascal; % pressure of active gas (O2 or H2)
            cOH  = 1000*mol/(meter^3); 
            T    = 333.15*Kelvin;
            R    = con.R;
            F    = con.F;
            
            state.(inm) = model.(inm).setupOHconcentration();
            
            nc = model.G.cells.num;
            state.T = T*ones(nc, 1);
            
            state = model.dispatchTemperature(state);
            
            liqrho = PorousTransportLayer.density(cOH, T);
            
            eldes = {oer, her};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};
                nc = model.(elde).(ptl).G.cells.num;

                fun = @(s) leverett(model.(elde).(ptl).leverettCoefficients, s); % Define Leverett function handle
                sLiquid = fzero(fun, 0.7); % Solve equilibrium liquid saturation

                svf = model.(elde).(ptl).solidVolumeFraction;
                lvf = sLiquid.*(1 - svf);
                
                state.(elde).(ptl).liqeps    = sLiquid.*(1 - svf);
                state.(elde).(ptl).OHceps    = cOH.*lvf;
                state.(elde).(ptl).liqrhoeps = liqrho*lvf;

                state.(elde).(ptl).phi = zeros(nc, 1);
                
                state.(elde).(ptl).liqrho = liqrho;
                state.(elde).(ptl) = model.(elde).(ptl).updateVolumeFractions(state.(elde).(ptl));
                state.(elde).(ptl) = model.(elde).(ptl).updateOHconcentration(state.(elde).(ptl));
                state.(elde).(ptl) = model.(elde).(ptl).updateMolality(state.(elde).(ptl));
                state = model.dispatchTemperature(state);
                state.(elde) = model.(elde).dispatchTemperature(state.(elde));
                state.(elde).(ptl) = model.(elde).(ptl).updateVaporPressure(state.(elde).(ptl));

                H2Ovp = state.(elde).(ptl).vaporPressure;
                H2Ogrho = H2Ovp*model.(elde).(ptl).sp.H2O.MW./(con.R*T);
                gvf = (1 - lvf - model.(elde).(ptl).solidVolumeFraction); % Gas volume fraction
                
                state.(elde).(ptl).H2Ogasrhoeps = H2Ogrho.*gvf;

            end

            switch elde
                
              case her
                
                H2p = pGas;
                H2rho = H2p.*model.(elde).(ptl).sp.H2.MW / (con.R * T);
                state.(elde).(ptl).H2rhoeps = H2rho*gvf;
                state.(elde).(ctl).E = 0;

              case oer

                O2p = pGas;
                O2rho = O2p.*model.(elde).(ptl).sp.O2.MW / (con.R * T);
                state.(elde).(ptl).O2rhoeps = O2rho*gvf;

                state.(elde).(ptl) = model.(elde).(ptl).updateVolumeFractions(state.(elde).(ptl));
                state.(elde).(ptl) = model.(elde).(ptl).updateOHconcentration(state.(elde).(ptl));
                state.(elde).(ptl) = model.(elde).(ptl).updateLiquidDensity(state.(elde).(ptl));
                state.(elde).(ptl) = model.(elde).(ptl).updateMolality(state.(elde).(ptl));
                state              = model.dispatchTemperature(state);
                state.(elde)       = model.(elde).dispatchTemperature(state.(elde));
                state.(elde).(ptl) = model.(elde).(ptl).updateWaterActivity(state.(elde).(ptl));
                state.(elde).(ptl) = model.(elde).(ptl).updateGasPressure(state.(elde).(ptl));
                state.(elde)       = model.(elde).dispatchToCatalystAndExchangeLayers(state.(elde));
                state.(elde).(ctl) = model.(elde).(ctl).updateEelyte(state.(elde).(ctl));
                
                Eelyte = state.(elde).(ctl).Eelyte;
                
                state.(elde).(ctl).E = Eelyte(1);
                
              otherwise
                
                error('electrode not recognized');
                
            end



            nc = model.(inm).G.cells.num;

            % We use water activity in oer. It has been already computed above.

            aw = state.(oer).(ptl).H2Oa(1);
            cH2O = IonomerMembrane.groupHydration(model.(inm), aw, T);
            model.(inm).H2O.c0 = cH2O/aw;
            state.(inm).H2Oceps = cH2O.*model.(inm).volumeFraction;

            state.(inm) = model.(inm).setupOHconcentration(state.(inm));
            cT  = state.(inm).cOH(1);
            phi = - R*T/F*log(cOH/cT);

            state.(inm).phi = phi*ones(nc, 1);

            bd = 'Boundary';
            
            state = model.evalVarNames(state, '.*compGasMasses');
            state = model.evalVarNames(state, 'Por.*phasePressures');
            state = model.evalVarNames(state, 'Layer.(?!Boundary).*liqrho$');
            
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                coupterm = model.(elde).(ptl).externalCouplingTerm;
                coupcells = coupterm.couplingcells;
                nc = numel(coupcells);

                gasInd = model.(elde).(ptl).gasInd;
                gvf = state.(elde).(ptl).volumeFractions{model.(elde).(ptl).phaseInd.gas};
                for igas = 1 :  gasInd.ngas
                    compgasdensity = state.(elde).(ptl).compGasMasses{igas}./gvf;
                    compgasdensity = compgasdensity(coupcells);
                    state.(elde).(ptl).(bd).gasDensities{igas} = compgasdensity;
                    controlValues.gasDensities{igas} = compgasdensity;
                end

                mobPhaseInd = model.(elde).(ptl).mobPhaseInd;
                for imobphase = 1 : mobPhaseInd.nmobphase
                    iphase = mobPhaseInd.phaseMap(imobphase);
                    p = state.(elde).(ptl).phasePressures{iphase};
                    p = p(coupcells);
                    state.(elde).(ptl).(bd).phasePressures{iphase} = p;
                    controlValues.mobilePhasePressures{imobphase} = p;
                end

                liqrho = state.(elde).(ptl).liqrho;
                liqrho = liqrho(coupcells);
                state.(elde).(ptl).(bd).liqrho = liqrho;
                controlValues.liqrho = liqrho;

                liquidInd = model.(elde).(ptl).liquidInd;
                cOH = state.(elde).(ptl).concentrations{liquidInd.OH};
                cOH = cOH(coupcells);
                state.(elde).(ptl).(bd).cOH = cOH;
                controlValues.cOH = cOH;

                model.(elde).(ptl).(bd).controlValues = controlValues;
                
            end
            
        end

        function state = setupControl(model, state, drivingForces)
            
            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';
            ctl = 'CatalystLayer';

            time = state.time;
            
            controlI = drivingForces.src(time);

            controlEqs{1} = state.(oer).(ctl).I - controlI;
            controlEqs{2} = state.(her).(ctl).E; % we impose zero potential at cathode

            state.controlEqs = controlEqs;
            
        end
        
        function state = dispatchTemperature(model, state)
            
            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';
            inm = 'IonomerMembrane';
            
            state.(oer).T = state.T(model.(oer).G.mappings.cellmap);
            state.(her).T = state.T(model.(her).G.mappings.cellmap);
            state.(inm).T = state.T(model.(inm).G.mappings.cellmap);
            
        end
        
        function state = dispatchIonomerToReactionLayers(model, state)

            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';
            inm = 'IonomerMembrane';

            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';

            eldes    = {her, oer};
            layers   = {ctl, exl};
            varnames = {'H2OaInmr', 'cOHinmr', 'phiInmr'};

            for ielde = 1 : numel(eldes)
                
                elde = eldes{ielde};

                % setup initial value for the variable that also takes care of AD.
                nc = model.(elde).(ctl).G.cells.num;
                initval = nan(nc, 1);
                [adsample, isAD] = getSampleAD(state.(inm).phi);
                if isAD
                    initval = adsample.convertDouble(initval);
                end

                coupterms = model.couplingTerms;
                coupnames = model.couplingNames;
                coupname  = sprintf('%s-%s', elde, inm);
                coupterm  = getCoupTerm(coupterms, coupname, coupnames);
                coupcells = coupterm.couplingcells;

                for ilayer = 1 : numel(layers)
                    layer = layers{ilayer};
                    for ivar = 1 : numel(varnames)
                        varname = varnames{ivar};
                        state.(elde).(layer).(varname) = initval;
                    end
                    state.(elde).(layer).H2OaInmr(coupcells(:, 1)) = state.(inm).H2Oa(coupcells(:, 2));
                    state.(elde).(layer).cOHinmr(coupcells(:, 1))  = state.(inm).cOH(coupcells(:, 2));
                    state.(elde).(layer).phiInmr(coupcells(:, 1))  = state.(inm).phi(coupcells(:, 2));
                end
                
            end
            
        end

        function state = updateIonomerSources(model, state)
            
            inm = 'IonomerMembrane';
            her = 'HydrogenEvolutionElectrode';
            oer = 'OxygenEvolutionElectrode';

            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';

            % initialize sources (done in a way that takes care of AD, meaning that OHsource and H2OSource inherits AD
            % structure from state.(inm).H2Oceps)
            OHsource  = 0*state.(inm).H2Oceps;
            H2OSource = 0*state.(inm).H2Oceps;

            eldes = {her, oer};

            for ielde = 1 : numel(eldes)
                
                elde = eldes{ielde};

                coupterms = model.couplingTerms;
                coupnames = model.couplingNames;
                vols      = model.G.cells.volumes;

                coupname  = sprintf('%s-%s', elde, inm);
                coupterm  = getCoupTerm(coupterms, coupname, coupnames);
                coupcells = coupterm.couplingcells;
                vols      = vols(coupcells(:, 2));
                
                inmrOHsource = state.(elde).(ctl).inmrReactionRate(coupcells(:, 1));
                OHexchR      = state.(elde).(exl).OHexchangeRate(coupcells(:, 1));
                H2OexchR     = state.(elde).(exl).H2OexchangeRate(coupcells(:, 1));
                
                OHsource(coupcells(:, 2))  = vols.*(inmrOHsource - OHexchR);
                H2OSource(coupcells(:, 2)) = - vols.*H2OexchR;
                
            end
            
            state.(inm).OHsource  = OHsource;
            state.(inm).H2OSource = H2OSource;
            
        end

        function primaryvarnames = getPrimaryVariables(model)
            
            primaryvarnames = model.primaryVarNames;
            
        end


        function cleanState = addStaticVariables(model, cleanState, state)

            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';
            ptl = 'PorousTransportLayer';
            inm = 'IonomerMembrane';
            
            cleanState = addStaticVariables@BaseModel(model, cleanState, state);
            cleanState.T = state.T;            
            cleanState.(inm).jBcSource = 0;
            cleanState.(oer).(ptl).jBcSource = 0;
            cleanState.(her).(ptl).jBcSource = 0;
            
        end

        
        function forces = getValidDrivingForces(model)
            
            forces = getValidDrivingForces@PhysicalModel(model);
            forces.src = [];
            
        end

        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
            opts = struct('ResOnly', false, 'iteration', 0); 
            opts = merge_options(opts, varargin{:});
            
            inm = 'IonomerMembrane';
            her = 'HydrogenEvolutionElectrode';
            oer = 'OxygenEvolutionElectrode';
            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';
            ptl = 'PorousTransportLayer';
            bd  = 'Boundary';

            time = state0.time + dt;
            if(not(opts.ResOnly))
                state = model.initStateAD(state);
            end

            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end
            
            %% Set up the governing equations
            
            eqs = {};
            names = {};
            
            eqs{end + 1}   = state.(inm).chargeCons;
            names{end + 1} = 'inm_chargeCons';
            eqs{end + 1}   = state.(inm).H2OmassCons;
            names{end + 1} = 'inm_H2OmassCons';
            eqs{end + 1}   = state.(oer).(ptl).(bd).bcEquations{1};
            names{end + 1} = 'oer_ptl_bd_bcEquations_1';
            eqs{end + 1}   = state.(oer).(ptl).(bd).bcEquations{2};
            names{end + 1} = 'oer_ptl_bd_bcEquations_2';
            eqs{end + 1}   = state.(oer).(ptl).(bd).bcEquations{3};
            names{end + 1} = 'oer_ptl_bd_bcEquations_3';
            eqs{end + 1}   = state.(oer).(ptl).(bd).bcEquations{4};
            names{end + 1} = 'oer_ptl_bd_bcEquations_4';
            eqs{end + 1}   = state.(oer).(ptl).(bd).bcControlEquations{1};
            names{end + 1} = 'oer_ptl_bd_bcControlEquations_1';
            eqs{end + 1}   = state.(oer).(ptl).(bd).bcControlEquations{2};
            names{end + 1} = 'oer_ptl_bd_bcControlEquations_2';
            eqs{end + 1}   = state.(oer).(ptl).chargeCons;
            names{end + 1} = 'oer_ptl_chargeCons';
            eqs{end + 1}   = state.(oer).(ptl).compGasMassCons{1};
            names{end + 1} = 'oer_ptl_compGasMassCons_1';
            eqs{end + 1}   = state.(oer).(ptl).compGasMassCons{2};
            names{end + 1} = 'oer_ptl_compGasMassCons_2';
            eqs{end + 1}   = state.(oer).(ptl).liquidMassCons;
            names{end + 1} = 'oer_ptl_liquidMassCons';
            eqs{end + 1}   = state.(oer).(ptl).OHMassCons;
            names{end + 1} = 'oer_ptl_OHMassCons';
            eqs{end + 1}   = state.(oer).(ptl).liquidStateEquation;
            names{end + 1} = 'oer_ptl_liquidStateEquation';
            eqs{end + 1}   = state.(her).(ptl).(bd).bcEquations{1};
            names{end + 1} = 'her_ptl_bd_bcEquations_1';
            eqs{end + 1}   = state.(her).(ptl).(bd).bcEquations{2};
            names{end + 1} = 'her_ptl_bd_bcEquations_2';
            eqs{end + 1}   = state.(her).(ptl).(bd).bcEquations{3};
            names{end + 1} = 'her_ptl_bd_bcEquations_3';
            eqs{end + 1}   = state.(her).(ptl).(bd).bcEquations{4};
            names{end + 1} = 'her_ptl_bd_bcEquations_4';
            eqs{end + 1}   = state.(her).(ptl).(bd).bcControlEquations{1};
            names{end + 1} = 'her_ptl_bd_bcControlEquations_1';
            eqs{end + 1}   = state.(her).(ptl).(bd).bcControlEquations{2};
            names{end + 1} = 'her_ptl_bd_bcControlEquations_2';
            eqs{end + 1}   = state.(her).(ptl).chargeCons;
            names{end + 1} = 'her_ptl_chargeCons';
            eqs{end + 1}   = state.(her).(ptl).compGasMassCons{1};
            names{end + 1} = 'her_ptl_compGasMassCons_1';
            eqs{end + 1}   = state.(her).(ptl).compGasMassCons{2};
            names{end + 1} = 'her_ptl_compGasMassCons_2';
            eqs{end + 1}   = state.(her).(ptl).liquidMassCons;
            names{end + 1} = 'her_ptl_liquidMassCons';
            eqs{end + 1}   = state.(her).(ptl).OHMassCons;
            names{end + 1} = 'her_ptl_OHMassCons';
            eqs{end + 1}   = state.(her).(ptl).liquidStateEquation;
            names{end + 1} = 'her_ptl_liquidStateEquation';
            eqs{end + 1}   = state.controlEqs{1};
            names{end + 1} = 'controlEqs_1';
            eqs{end + 1}   = state.controlEqs{2};
            names{end + 1} = 'controlEqs_2';

            neq = numel(eqs);
            
            types = repmat({'cell'}, 1, neq);
            
            primaryVars = model.getPrimaryVariables();

            %% setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
        
        end
        
    end
    
end


