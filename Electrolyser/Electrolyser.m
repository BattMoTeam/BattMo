classdef Electrolyser < BaseModel
    
    properties
        
        con = PhysicalConstants();

        % Components
        IonomerMembrane
        HydrogenEvolutionElectrode
        OxygenEvolutionElectrode
        
        couplingTerms
        couplingNames
        
    end
    
    methods
        
        function model = Electrolyser(paramobj)

            model = model@BaseModel();

            model.HydrogenEvolutionElectrode = EvolutionElectrode(paramobj.HydrogenEvolutionElectrode);
            model.OxygenEvolutionElectrode   = EvolutionElectrode(paramobj.OxygenEvolutionElectrode);
            model.IonomerMembrane            = IonomerMembrane(paramobj.IonomerMembrane);
            
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {'T'};

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
                inputvarnames{end + 1} = {elde, ctl, 'elyteReactionRate'};
                inputvarnames{end + 1} = {elde, ctl, 'inmrReactionRate'};
                inputvarnames{end + 1} = {elde, exl, 'OHexchangeRate'};
                inputvarnames{end + 1} = {elde, exl, 'H2OexchangeRate'};
            end
            model = model.registerPropFunction({{inm, 'H2OSource'}, fn, inputvarnames});
            model = model.registerPropFunction({{inm, 'OHSource'}, fn, inputvarnames});


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
            
            model = model.registerStaticVarName('T');
            
        end

        
        function initstate = setupInitialState(model)

        end

        function state = dispatchTemperature(model, state)
            
            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';
            inm = 'IonomerMembrane';
            
            state.(oer).T = state.T(model.(oer).G.cellmap);
            state.(her).T = state.T(model.(her).G.cellmap);
            state.(inm).T = state.T(model.(inm).G.cellmap);
            
        end
        
        function state = dispatchIonomerToReactionLayers(model, state)

            inm = 'IonomerMembrane';
            her = 'HydrogenEvolutionElectrode';
            oer = 'OxygenEvolutionElectrode';

            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';

            eldes = {her, oer};

            for ielde = 1 : numel(eldes)
                
                elde = eldes{ielde};
                coupterms = model.couplingTerms;
                coupnames = model.couplingNames;
                coupterm  = getCoupTerm(coupterms, {inm, elde}, coupnames);
                coupcells = coupterm.couplingcells;

                for ilayer = 1 : numel(layers)
                    layer = layers{ilayer};
                    state.(elde).(layer).H2OaInmr(coupcells(:, 2)) = state.(inm).H2Oa(coupcells(:, 1));
                    state.(elde).(layer).cOHinmr(coupcells(:, 2))  = state.(inm).cOH(coupcells(:, 1));
                    state.(elde).(layer).phiInmr(coupcells(:, 2))  = state.(inm).phi(coupcells(:, 1));
                end
                
            end
            
        end

        function state = updateIonomerSources(model, state)
            
            inm = 'IonomerMembrane';
            her = 'HydrogenEvolutionElectrode';
            oer = 'OxygenEvolutionElectrode';

            ctl = 'CatalystLayer';
            exl = 'ExchangeLayer';

            % initialize sources (done in a way that takes care of AD, meaning that OHSource and H2OSource are AD
            % variables is state.(inm).phi is)
            OHSource = 0*state.(inm).phi;
            H2OSource = 0*state.(inm).phi;

            eldes = {her, oer};

            for ielde = 1 : numel(eldes)
                
                elde = eldes{ielde};
                coupterms = model.couplingTerms;
                coupnames = model.couplingNames;
                coupterm  = getCoupTerm(coupterms, {inm, elde}, coupnames);
                coupcells = coupterm.couplingcells;

                reacR    = state.(elde).(ctl).inmrReactionRate(coupcells(:, 2));
                OHexchR  = state.(elde).(exl).OHexchangeRate(coupcells(:, 2));
                H2OexchR = state.(elde).(exl).H2OexchangeRate(coupcells(:, 2));
                
                leps  = model.(inm).liquidVolumeFraction(coupcells(:, 2));
                
                OHSource(coupcells(:, 1))  = -2.*reacR./(n.*F) - OHexchR.*leps;
                H2OSource(coupcells(:, 1)) = -reacR./(n.*F) - H2OexchR;
                
            end
            
            state.(inm).OHSource  = OHSource;
            state.(inm).H2OSource = H2OSource;
            
        end
        
        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
            opts = struct('ResOnly', false, 'iteration', 0); 
            opts = merge_options(opts, varargin{:});
            
            time = state0.time + dt;
            if(not(opts.ResOnly))
                state = model.initStateAD(state);
            end
            
            %% for now temperature and SOC are kept constant
            nc = model.G.cells.num;
            state.SOC = model.SOC*ones(nc, 1);
            
            % Shortcuts used in this function
            battery = model;
            ne      = 'NegativePorousTransportLayer';
            pe      = 'PositivePorousTransportLayer';
            eac     = 'PorousTransportLayerActiveComponent';
            cc      = 'CurrentCollector';
            elyte   = 'Electrolyte';
            am      = 'ActiveMaterial';
            thermal = 'ThermalModel';
            
            electrodes = {ne, pe};
            electrodecomponents = {eac, cc};

            %% Synchronization across components

            % temperature
            state = battery.updateTemperature(state);
            
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                % potential and concentration between active material and electode active component
                state.(elde).(eac) = battery.(elde).(eac).updatePhi(state.(elde).(eac));
                state.(elde).(eac) = battery.(elde).(eac).updateChargeCarrier(state.(elde).(eac));
                
            end
            
            %% Accumulation terms

            state = battery.updateAccumTerms(state, state0, dt);

            %% Update Electrolyte -> PorousTransportLayers coupling 
            
            state = battery.updatePorousTransportLayerCoupling(state); 

            %% Update reaction rates in both electrodes

            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateDiffusionConductivityCoefficients(state.(elde).(eac).(am));
                state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateOCP(state.(elde).(eac).(am));
                state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateReactionRate(state.(elde).(eac).(am));
            end

            %% Update PorousTransportLayers -> Electrolyte  coupling

            state = battery.updateElectrolyteCoupling(state);
            
            %% Update coupling within electrodes and external coupling
            
            state.(ne) = battery.(ne).updateCoupling(state.(ne));
            state.(pe) = battery.(pe).updateCoupling(state.(pe));

            state.(ne).(eac) = battery.(ne).(eac).updatejBcSource(state.(ne).(eac));
            state.(pe).(eac) = battery.(pe).(eac).updatejBcSource(state.(pe).(eac));
            
            state = model.setupExternalCouplingNegativePorousTransportLayer(state);
            state = model.setupExternalCouplingPositivePorousTransportLayer(state);
            
            state.(ne).(cc) = battery.(ne).(cc).updatejBcSource(state.(ne).(cc));
            state.(pe).(cc) = battery.(pe).(cc).updatejBcSource(state.(pe).(cc));
                        
            %% elyte charge conservation

            state.(elyte) = battery.(elyte).updateCurrentBcSource(state.(elyte));
            state.(elyte) = battery.(elyte).updateChemicalCurrent(state.(elyte));
            state.(elyte) = battery.(elyte).updateCurrent(state.(elyte));
            state.(elyte) = battery.(elyte).updateChargeConservation(state.(elyte));

            %% PorousTransportLayers charge conservation - Active material part

            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                state.(elde).(eac) = battery.(elde).(eac).updateIonAndCurrentSource(state.(elde).(eac));
                state.(elde).(eac) = battery.(elde).(eac).updateCurrent(state.(elde).(eac));
                state.(elde).(eac) = battery.(elde).(eac).updateChargeConservation(state.(elde).(eac));
            end
            
            %% elyte mass conservation

            state.(elyte) = battery.(elyte).updateDiffusionCoefficient(state.(elyte));
            state.(elyte) = battery.(elyte).updateChargeCarrierFlux(state.(elyte));
            state.(elyte) = battery.(elyte).updateMassConservation(state.(elyte));

            
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                
                %% PorousTransportLayers mass conservation
                state.(elde).(eac) = battery.(elde).(eac).updateDiffusionCoefficient(state.(elde).(eac));
                state.(elde).(eac) = battery.(elde).(eac).updateChargeCarrierFlux(state.(elde).(eac));
                state.(elde).(eac) = battery.(elde).(eac).updateMassConservation(state.(elde).(eac));
                
                %% PorousTransportLayers charge conservation - current collector part
                state.(elde).(cc) = battery.(elde).(cc).updateCurrent(state.(elde).(cc));
                state.(elde).(cc) = battery.(elde).(cc).updateChargeConservation(state.(elde).(cc));

            end
            
            %% update Thermal source term from electrical resistance

            state = battery.updateThermalOhmicSourceTerms(state);
            state = battery.updateThermalChemicalSourceTerms(state);
            state = battery.updateThermalReactionSourceTerms(state);
           
            state.(thermal) = battery.(thermal).updateHeatSourceTerm(state.(thermal));
            state.(thermal) = battery.(thermal).updateThermalBoundarySourceTerms(state.(thermal));
            
            %% update Accumulation terms for the energy equation
            
            state = battery.updateThermalAccumTerms(state, state0, dt);
            
            %% Update energy conservation residual term
            
            state.(thermal) = model.(thermal).updateEnergyConservation(state.(thermal));
            
            
            %% Set up the governing equations
            
            eqs = {};
            
            %% We collect mass and charge conservation equations for the electrolyte and the electrodes

            eqs{end + 1} = state.(elyte).massCons;
            eqs{end + 1} = state.(elyte).chargeCons;
            
            eqs{end + 1} = state.(ne).(eac).massCons;
            eqs{end + 1} = state.(ne).(eac).chargeCons;
            
            eqs{end + 1} = state.(pe).(eac).massCons;
            eqs{end + 1} = state.(pe).(eac).chargeCons;
            
            eqs{end + 1} = state.(ne).(cc).chargeCons;
            eqs{end + 1} = state.(pe).(cc).chargeCons;
            
            eqs{end + 1} = state.(thermal).energyCons;
            
            %% We setup and add the control equation (fixed total current at PositiveCurrentCollector)
            
            src = drivingForces.src(time);
            coupterm = battery.(pe).(cc).couplingTerm;
            faces = coupterm.couplingfaces;
            bcval = state.(pe).(cc).E;
            cond_pcc = battery.(pe).(cc).EffectiveElectricalConductivity;
            [trans_pcc, cells] = battery.(pe).(cc).operators.harmFaceBC(cond_pcc, faces);
            control = sum(trans_pcc.*(state.(pe).(cc).phi(cells) - bcval)) - src;
            
            eqs{end + 1} = - control;

            %% Give type and names to equations and names of the primary variables (for book-keeping)
            
            types = {'cell','cell','cell','cell',...
                     'cell','cell','cell','cell','cell','cell'};
            
            names = {'elyte_massCons'   , ...
                     'elyte_chargeCons' , ...
                     'ne_eac_massCons'  , ...
                     'ne_eac_chargeCons', ...
                     'pe_eac_massCons'  , ...
                     'pe_eac_chargeCons', ...
                     'ne_cc_chargeCons' , ...
                     'pe_cc_chargeCons' , ...
                     'energyCons'       , ...
                     'control'};
            
            primaryVars = model.getPrimaryVariables();

            %% setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
        
        end

        function state = updateTemperature(model, state)
            
            elyte = 'Electrolyte';
            ne    = 'NegativePorousTransportLayer';
            pe    = 'PositivePorousTransportLayer';
            eac   = 'PorousTransportLayerActiveComponent';
            cc    = 'CurrentCollector';
            thermal = 'ThermalModel';
            
            % (here we assume that the ThermalModel has the "parent" grid)
            state.(elyte).T    = state.(thermal).T(model.(elyte).G.mappings.cellmap);
            state.(ne).(eac).T = state.(thermal).T(model.(ne).(eac).G.mappings.cellmap);
            state.(ne).(cc).T  = state.(thermal).T(model.(ne).(cc).G.mappings.cellmap);
            state.(pe).(eac).T = state.(thermal).T(model.(pe).(eac).G.mappings.cellmap);
            state.(pe).(cc).T  = state.(thermal).T(model.(pe).(cc).G.mappings.cellmap);
            
            % Update temperature in the active materials of the electrodes.
            state.(ne).(eac) = model.(ne).(eac).updateTemperature(state.(ne).(eac));
            state.(pe).(eac) = model.(pe).(eac).updateTemperature(state.(pe).(eac));
            
        end
        
        
        function state = updateElectrolyteCoupling(model, state)
        % Setup the electrolyte coupling by adding ion sources from the electrodes
        % shortcuts:
        % c_source : Source term for charge carrier.
                        
            battery = model;
            elyte = 'Electrolyte';
            ne    = 'NegativePorousTransportLayer';
            pe    = 'PositivePorousTransportLayer';
            am    = 'ActiveMaterial';
            eac   = 'PorousTransportLayerActiveComponent';
            
            vols = battery.(elyte).G.cells.volumes;
            F = battery.con.F;
            
            ccSourceName = battery.(elyte).chargeCarrierSourceName;
            couplingterms = battery.couplingTerms;

            elyte_c_source = zeros(battery.(elyte).G.cells.num, 1);
            elyte_e_source = zeros(battery.(elyte).G.cells.num, 1);
            
            % setup AD 
            phi = state.(elyte).phi;
            if isa(phi, 'ADI')
                adsample = getSampleAD(phi);
                adbackend = model.AutoDiffBackend;
                elyte_c_source = adbackend.convertToAD(elyte_c_source, adsample);
            end
            
            coupnames = model.couplingNames;
            
            ne_R = state.(ne).(eac).(am).R;
            coupterm = getCoupTerm(couplingterms, 'NegativePorousTransportLayer-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = ne_R.*vols(elytecells); % we divide with F later
            
            pe_R = state.(pe).(eac).(am).R;
            coupterm = getCoupTerm(couplingterms, 'PositivePorousTransportLayer-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = pe_R.*vols(elytecells);
            
            elyte_e_source = elyte_c_source.*battery.(elyte).sp.z(1); % we divide with F later
            
            state.Electrolyte.(ccSourceName) = elyte_c_source/F; 
            state.Electrolyte.eSource = elyte_e_source;
            
        end
        
        function state = updateAccumTerms(model, state, state0, dt)
                    
            elyte = 'Electrolyte';
            ne    = 'NegativePorousTransportLayer';
            pe    = 'PositivePorousTransportLayer';
            am    = 'ActiveMaterial';
            eac   = 'PorousTransportLayerActiveComponent';
            
            ccAccumName = model.(elyte).chargeCarrierAccumName;
            
            cdotcc  = (state.(elyte).cs{1} - state0.(elyte).cs{1})/dt;
            effectiveVolumes = model.(elyte).volumeFraction.*model.(elyte).G.cells.volumes;
            ccAccum  = effectiveVolumes.*cdotcc;
            state.(elyte).(ccAccumName) = ccAccum;
            
            names = {ne, pe};
            for i = 1 : numel(names)
                elde = names{i}; % electrode name
                cdotcc   = (state.(elde).(eac).c - state0.(elde).(eac).c)/dt;
                effectiveVolumes = model.(elde).(eac).volumeFraction.*model.(elde).(eac).G.cells.volumes;
                ccAccum  = effectiveVolumes.*cdotcc;
                state.(elde).(eac).(ccAccumName) = ccAccum;
            end
            
        end

        
        function state = updateThermalAccumTerms(model, state, state0, dt)
                    
            thermal = 'ThermalModel';
         
            hcap = model.(thermal).EffectiveHeatCapacity;
            
            T = state.(thermal).T;
            T0 = state0.(thermal).T;

            % (here we assume that the ThermalModel has the "parent" grid)
            vols = model.G.cells.volumes;
            
            state.(thermal).accumHeat = hcap.*vols.*(T - T0)/dt;
            
        end


        function state = updateThermalOhmicSourceTerms(model, state)
        % reference Latz et al (ref1 in reference list)

            ne      = 'NegativePorousTransportLayer';
            pe      = 'PositivePorousTransportLayer';
            eac     = 'PorousTransportLayerActiveComponent';
            cc      = 'CurrentCollector';
            elyte   = 'Electrolyte';
            thermal = 'ThermalModel';
            
            eldes = {ne, pe}; % electrodes
            
            nc = model.G.cells.num;
            
            src = zeros(nc, 1);
            
            T = state.(thermal).T;
            if isa(T, 'ADI')
                adsample = getSampleAD(T);
                adbackend = model.AutoDiffBackend;
                src = adbackend.convertToAD(src, adsample);
            end
            
            for ind = 1 : numel(eldes)

                elde = eldes{ind};
                
                cc_model = model.(elde).(cc);
                cc_map   = cc_model.G.mappings.cellmap;
                cc_j     = state.(elde).(cc).j;
                cc_econd = cc_model.EffectiveElectricalConductivity;
                
                cc_src = computeCellFluxNorm(cc_model, cc_j, 1, cc_econd); % note that we send volumeFraction = 1
                src(cc_map) = src(cc_map) + cc_src;

                eac_model = model.(elde).(eac);
                eac_map   = eac_model.G.mappings.cellmap;
                eac_j     = state.(elde).(eac).j;
                eac_econd = eac_model.EffectiveElectricalConductivity;
                eac_vf    = eac_model.volumeFraction;
                
                eac_src = computeCellFluxNorm(eac_model, eac_j, eac_vf, eac_econd);
                src(eac_map) = src(eac_map) + eac_src;
                
            end

            % Electrolyte
            elyte_model = model.(elyte);
            elyte_map   = elyte_model.G.mappings.cellmap;
            elyte_j     = state.(elyte).j;
            elyte_econd = state.(elyte).conductivity;
            elyte_vf    = elyte_model.volumeFraction;            
            
            elyte_src = computeCellFluxNorm(elyte_model, elyte_j, elyte_vf, elyte_econd);
            src(elyte_map) = src(elyte_map) + elyte_src;
            
            state.(thermal).jHeatOhmSource = src;
            
            % For simplicity, we set the boundary heat transfer term equal to zero here.
            state.(thermal).jHeatBcSource = 0;
            
        end
        
        function state = updateThermalChemicalSourceTerms(model, state)
        % reference Latz et al (ref1 in reference list)            
            
            elyte = 'Electrolyte';
            thermal = 'ThermalModel';
              
            % prepare term
            nc = model.G.cells.num;
            src = zeros(nc, 1);
            T = state.(thermal).T;
            if isa(T, 'ADI')
                adsample = getSampleAD(T);
                adbackend = model.AutoDiffBackend;
                src = adbackend.convertToAD(src, adsample);
            end            

            % Compute chemical heat source in electrolyte
            dmudcs = state.(elyte).dmudcs; % Derivative of chemical potential with respect to concentration
            D = state.(elyte).D;           % Diffusion coefficient 
            Ddc = state.(elyte).diffFlux;  % Diffusion flux (-D*grad(c))
            
            elyte_map = model.(elyte).G.mappings.cellmap;
            vf = model.(elyte).volumeFraction;
            
            elyte_src = computeCellFluxNorm(model.(elyte), Ddc, vf, D);
            % This is a bit hacky for the moment (we should any way consider all the species)
            elyte_src = dmudcs{1}.*elyte_src;
            
            % map to source term at battery level
            src(elyte_map) = src(elyte_map) + elyte_src;
            
            state.(thermal).jHeatChemicalSource = src;
            
        end
        
        
        function state = updateThermalReactionSourceTerms(model, state)
        % reference Latz et al (ref1 in reference list)            
            
            ne      = 'NegativePorousTransportLayer';
            pe      = 'PositivePorousTransportLayer';
            eac     = 'PorousTransportLayerActiveComponent';
            am      = 'ActiveMaterial';
            thermal = 'ThermalModel';
            
            eldes = {ne, pe}; % electrodes
            
            nc = model.G.cells.num;
            
            src = zeros(nc, 1);
            
            T = state.(thermal).T;
            if isa(T, 'ADI')
                adsample = getSampleAD(T);
                adbackend = model.AutoDiffBackend;
                src = adbackend.convertToAD(src, adsample);
            end
            
            for ind = 1 : numel(eldes)

                elde = eldes{ind};
                am_model = model.(elde).(eac).(am);
                am_map   = am_model.G.mappings.cellmap;
                R = state.(elde).(eac).(am).R;
                eta = state.(elde).(eac).(am).eta;
                vols = model.(elde).(eac).G.cells.volumes;
                am_src = vols.*R.*eta;
                
                src(am_map) = src(am_map) + am_src;
                
            end

            % We multiply by volumes
            src = model.G.cells.volumes.*src;
            
            state.(thermal).jHeatReactionSource = src;

        end
        
        
        function state = updatePorousTransportLayerCoupling(model, state)
        % Setup electrode coupling by updating the potential and concentration of the electrolyte in the active component of the
        % electrode. There, those quantities are considered as input and used to compute the reaction rate.
        %
        %
        % WARNING : at the moment, we do not pass the concentrations
        %
        % shortcuts:
        % elyte : Electrolyte
        % neac  : NegativePorousTransportLayer.PorousTransportLayerActiveComponent 
        % peac  : PositivePorousTransportLayer.PorousTransportLayerActiveComponent
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativePorousTransportLayer';
            pe    = 'PositivePorousTransportLayer';
            am    = 'ActiveMaterial';
            eac   = 'PorousTransportLayerActiveComponent';
            cc    = 'CurrentCollector';
            
            eldes = {ne, pe};
            phi_elyte = state.(elyte).phi;
            c_elyte = state.(elyte).cs{1};
            
            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(bat.(elyte).G.mappings.cellmap) = (1 : bat.(elyte).G.cells.num)';

            for ind = 1 : numel(eldes)
                elde = eldes{ind};
                state.(elde).(eac).(am).phiElectrolyte = phi_elyte(elyte_cells(bat.(elde).(eac).G.mappings.cellmap));
                state.(elde).(eac).(am).cElectrolyte = c_elyte(elyte_cells(bat.(elde).(eac).G.mappings.cellmap));
            end
            
        end

        function state = setupExternalCouplingNegativePorousTransportLayer(model, state)
            
            ne = 'NegativePorousTransportLayer';
            cc = 'CurrentCollector';
           
            phi = state.(ne).(cc).phi;

            jExternal = model.(ne).(cc).setupExternalCoupling(phi, 0);
            
            state.(ne).(cc).jExternal = jExternal;
            
        end
        
        function state = setupExternalCouplingPositivePorousTransportLayer(model, state)
            
            pe = 'PositivePorousTransportLayer';
            cc = 'CurrentCollector';
           
            phi = state.(pe).(cc).phi;
            E = state.(pe).(cc).E;
            
            jExternal = model.(pe).(cc).setupExternalCoupling(phi, E);
            
            state.(pe).(cc).jExternal = jExternal;
            
        end

        
        function state = initStateAD(model,state)
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativePorousTransportLayer';
            pe    = 'PositivePorousTransportLayer';
            am    = 'ActiveMaterial';
            eac   = 'PorousTransportLayerActiveComponent';
            cc    = 'CurrentCollector';
            thermal = 'ThermalModel';
            
            adbackend = model.AutoDiffBackend();
            useMex=false;
            if(isprop(adbackend,'useMex'))
               useMex = adbackend.useMex; 
            end
            opts=struct('types',[1,1,2,2,3,3,4,5,6,7],'useMex',useMex);
            [state.(elyte).cs{1}  , ...
             state.(elyte).phi    , ...   
             state.(ne).(eac).c   , ...   
             state.(ne).(eac).phi , ...   
             state.(pe).(eac).c   , ...    
             state.(pe).(eac).phi , ...   
             state.(ne).(cc).phi  , ...    
             state.(pe).(cc).phi  , ...    
             state.(thermal).T    , ...
             state.(pe).(cc).E] = ...
                adbackend.initVariablesAD(...
                    state.(elyte).cs{1}  , ...
                    state.(elyte).phi    , ...   
                    state.(ne).(eac).c   , ...    
                    state.(ne).(eac).phi , ...   
                    state.(pe).(eac).c   , ...    
                    state.(pe).(eac).phi , ...   
                    state.(ne).(cc).phi  , ...    
                    state.(pe).(cc).phi  , ...    
                    state.(thermal).T    , ...
                    state.(pe).(cc).E, ...
                    opts); 
            % PRIMARY variables
        end
        
        function p = getPrimaryVariables(model)
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativePorousTransportLayer';
            pe    = 'PositivePorousTransportLayer';
            am    = 'ActiveMaterial';
            eac   = 'PorousTransportLayerActiveComponent';
            cc    = 'CurrentCollector';
            thermal = 'ThermalModel';
            
            p = {{elyte, 'cs', 1} , ...
                 {elyte, 'phi'}   , ...   
                 {ne, eac, 'c'}   , ...    
                 {ne, eac, 'phi'} , ...   
                 {pe, eac, 'c'}   , ...    
                 {pe, eac, 'phi'} , ...   
                 {ne, cc, 'phi'}  , ...    
                 {pe, cc, 'phi'}  , ...
                 {thermal, 'T'}   , ...
                 {pe, cc, 'E'}
                };
            
        end
        
        function state = setProp(model, state, names, val)
            if iscell(names) & (numel(names) > 1)
                name = names{1};
                names = names(2 : end);
                state.(name) = model.setProp(state.(name), names, val);
            elseif iscell(names) & (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    state{name} = val;
                else
                    state.(name) = val;
                end
            else
                error('format not recognized');                
            end                
        end
        
        function var = getProp(model, state, names)
            if iscell(names) && (numel(names) > 1)
                name = names{1};
                names = names(2 : end);
                var = model.getProp(state.(name), names);
            elseif iscell(names) & (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    var = state{name};
                else
                    var = state.(name);
                end
            else
                error('format not recognized');
            end                
        end
        
        function submod = getSubmodel(model, names)
            submod = model.(names{1});
            for i=2:numel(names)
                submod = submod.(names{i});
            end
        end


        
        function validforces = getValidDrivingForces(model)
            validforces=struct('src', [], 'stopFunction', []); 
        end
       function model = validateModel(model, varargin)
            mnames = {{'Electrolyte'}, ...
              {'PositivePorousTransportLayer','PorousTransportLayerActiveComponent'}, ...
              {'NegativePorousTransportLayer','PorousTransportLayerActiveComponent'}, ...
              {'NegativePorousTransportLayer','CurrentCollector'}, ...
              {'PositivePorousTransportLayer','CurrentCollector'}};
            model.Electrolyte.AutoDiffBackend=model.AutoDiffBackend;
            model.Electrolyte=model.Electrolyte.validateModel(varargin{:});
            model.PositivePorousTransportLayer.PorousTransportLayerActiveComponent.AutoDiffBackend= model.AutoDiffBackend;
            model.PositivePorousTransportLayer.PorousTransportLayerActiveComponent = model.PositivePorousTransportLayer.PorousTransportLayerActiveComponent.validateModel(varargin{:});
            model.NegativePorousTransportLayer.PorousTransportLayerActiveComponent.AutoDiffBackend= model.AutoDiffBackend;
            model.NegativePorousTransportLayer.PorousTransportLayerActiveComponent = model.NegativePorousTransportLayer.PorousTransportLayerActiveComponent.validateModel(varargin{:});
            model.NegativePorousTransportLayer.CurrentCollector.AutoDiffBackend= model.AutoDiffBackend;
            model.NegativePorousTransportLayer.CurrentCollector=model.NegativePorousTransportLayer.CurrentCollector.validateModel(varargin{:});
            model.PositivePorousTransportLayer.CurrentCollector.AutoDiffBackend=model.AutoDiffBackend;
            model.PositivePorousTransportLayer.CurrentCollector= model.PositiveElectrode.CurrentCollector.validateModel(varargin{:});
          %for i=1:numel(mnames)
          %    mname=mnames{i}
          %    submodel=model.getSubmodel(mname);
          %    submodel.AutoDiffBackend = model.AutoDiffBackend;
          %    submodel=submodel.validateModel(varargin{:});
          %    model  = model.setProp(model,mname,submodel);
          %end
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            p = model.getPrimaryVariables();
            for i=2:numel(dx)
                val = model.getProps(state, p{i});
                val = val + dx{i};
                state = model.setProp(state, p{i}, val);
            end
            %% not sure how to handle cells
            state.Electrolyte.cs{1} =  state.Electrolyte.cs{1} + dx{1};
            report = [];
        end
        
        
        function state = reduceState(model, state, removeContainers)
        % Reduce state to double (do not use property containers)
            state = value(state);
        end

        
    end
    
end


