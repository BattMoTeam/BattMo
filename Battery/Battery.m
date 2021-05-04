classdef Battery < PhysicalModel

    properties
        
        con = PhysicalConstants();

        % SOC
        SOC

        % Initial temperature
        initT

        % Input current
        J
        
        % Voltage cut
        Ucut

        % Components
        Electrolyte
        NegativeElectrode
        PositiveElectrode
        
        couplingTerms
        couplingNames
        
        mappings
        
    end
    
    methods
        
        function model = Battery(paramobj)
        % paramobj is instance of BatteryInputParams
        % Shorcuts used here
        % elyte : Electrolyte
        % ne : NegativeElectrode
        % pe : PositiveElectrode
            
            model = model@PhysicalModel([]);
            
            % OBS : All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            %% Setup the model using the input parameters
            fdnames = {'G'             , ...
                       'couplingTerms' , ...
                       'initT'         , ...
                       'SOC'           , ...
                       'J'             , ...
                       'Ucut'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
            % Assign the components : Electrolyte, NegativeElectrode, PositiveElectrode
            model.NegativeElectrode = model.setupElectrode(paramobj.ne);
            model.PositiveElectrode = model.setupElectrode(paramobj.pe);
            model.Electrolyte       = model.setupElectrolyte(paramobj.elyte);
            
            % setup couplingNames
            model.couplingNames = cellfun(@(x) x.name, model.couplingTerms, 'uniformoutput', false);
            
            % setup Electrolyte porosity
            model = model.setElectrolyteVolumeFraction();
            
            % setup Electrolyte effective thermal conductivity
            model = model.setElectrolyteEffectiveThermalConductivity();
            
            % setup mappings (electrodes -> electrolyte)
            model = model.setupMappings();
           
        end
        
        function electrode = setupElectrode(model, paramobj)
        % Standard setup (ActiveMaterial is specified in Electrode instantiations)
            electrode = Electrode(paramobj);
        end
        
        function electrode = setupElectrolyte(model, paramobj)
        % paramobj is instance of ElectrolyteInputParams
            switch paramobj.name
              case 'orgLiPF6'
                electrode = orgLiPF6(paramobj);
              otherwise
                error('electrolyte name not recognized');
            end
        end
        
        function model = setupMappings(model)
            
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            eac     = 'ElectrodeActiveComponent';
            cc      = 'CurrentCollector';
            elyte   = 'Electrolyte';
            
            
            G_elyte = model.(elyte).G;
            elytecelltbl.cells = (1 : G_elyte.cells.num)';
            elytecelltbl.globalcells = G_elyte.mappings.cellmap;
            elytecelltbl = IndexArray(elytecelltbl);

            eldes = {ne, pe};

            for ind = 1 : numel(eldes)

                elde = eldes{ind};
                G_elde  = model.(elde).(eac).G;
                clear eldecelltbl;
                eldecelltbl.cells = (1 : G_elde.cells.num)';
                eldecelltbl.globalcells = G_elde.mappings.cellmap;
                eldecelltbl = IndexArray(eldecelltbl);
                
                map = TensorMap();
                map.fromTbl = elytecelltbl;
                map.toTbl = eldecelltbl;
                map.replaceFromTblfds = {{'cells', 'elytecells'}};
                map.replaceToTblfds = {{'cells', 'eldecells'}};
                map.mergefds = {'globalcells'};
                
                mappings.(elde) = map.getDispatchInd();
                
            end
            
            model.mappings = mappings;
            
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
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            eac     = 'ElectrodeActiveComponent';
            cc      = 'CurrentCollector';
            elyte   = 'Electrolyte';
            am      = 'ActiveMaterial';
            
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

            %% Update Electrolyte -> Electrodes coupling 
            
            state = battery.updateElectrodeCoupling(state); 

            %% Update reaction rates in both electrodes

            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateDiffusionConductivityCoefficients(state.(elde).(eac).(am));
                state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateOCP(state.(elde).(eac).(am));
                state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateReactionRate(state.(elde).(eac).(am));
            end

            %% Update Electrodes -> Electrolyte  coupling

            state = battery.updateElectrolyteCoupling(state);
            
            %% Update coupling within electrodes and external coupling
            
            state.(ne) = battery.(ne).updateCoupling(state.(ne));
            state.(pe) = battery.(pe).updateCoupling(state.(pe));

            state.(ne).(eac) = battery.(ne).(eac).updatejBcSource(state.(ne).(eac));
            state.(pe).(eac) = battery.(pe).(eac).updatejBcSource(state.(pe).(eac));
            
            state = model.setupExternalCouplingNegativeElectrode(state);
            state = model.setupExternalCouplingPositiveElectrode(state);
            
            state.(ne).(cc) = battery.(ne).(cc).updatejBcSource(state.(ne).(cc));
            state.(pe).(cc) = battery.(pe).(cc).updatejBcSource(state.(pe).(cc));
                        
            %% elyte charge conservation

            state.(elyte) = battery.(elyte).updateCurrentBcSource(state.(elyte));
            state.(elyte) = battery.(elyte).updateChemicalCurrent(state.(elyte));
            state.(elyte) = battery.(elyte).updateCurrent(state.(elyte));
            state.(elyte) = battery.(elyte).updateChargeConservation(state.(elyte));

            %% Electrodes charge conservation - Active material part

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
                
                %% Electrodes mass conservation
                state.(elde).(eac) = battery.(elde).(eac).updateDiffusionCoefficient(state.(elde).(eac));
                state.(elde).(eac) = battery.(elde).(eac).updateChargeCarrierFlux(state.(elde).(eac));
                state.(elde).(eac) = battery.(elde).(eac).updateMassConservation(state.(elde).(eac));
                
                %% Electrodes charge conservation - current collector part
                state.(elde).(cc) = battery.(elde).(cc).updateCurrent(state.(elde).(cc));
                state.(elde).(cc) = battery.(elde).(cc).updateChargeConservation(state.(elde).(cc));

            end

            %% update Accumulation terms for the energy equation
            state = battery.updateEnergyAccumTerms(state, state0, dt);
            
            %% elyte energy conservation 
            state = battery.updateElyteHeatBcSource(state);
            state.(elyte) = battery.(elyte).updateOhmSource(state.(elyte));
            state.(elyte) = battery.(elyte).updateHeatSource(state.(elyte));
            state.(elyte) = battery.(elyte).updateEnergyConservation(state.(elyte));
            
            %% energy conservation in the electrodes
            for ind1 = 1 : numel(electrodes)
                elde = electrodes{ind1};
                state.(elde) = battery.(elde).updateTemperatureCoupling(state.(elde));
                for ind2 = 1 : numel(electrodecomponents)
                    comp = electrodecomponents{ind2};
                    state.(elde).(comp) = battery.(elde).(comp).updateOhmSource(state.(elde).(comp));
                    state.(elde).(comp) = battery.(elde).(comp).updateHeatSource(state.(elde).(comp));
                    state.(elde).(comp) = battery.(elde).(comp).updateEnergyConservation(state.(elde).(comp));
                end
                
            end

            %% energy conservation equation 
            state = battery.updateEnergyConservation(state);
            
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
            
            eqs{end + 1} = state.energyCons;
            
            %% We setup and add the control equation (fixed total current at PositiveCurrentCollector)
            
            src = drivingForces.src(time);
            coupterm = battery.(pe).(cc).couplingTerm;
            faces = coupterm.couplingfaces;
            bcval = state.(pe).(cc).E;
            cond_pcc = battery.(pe).(cc).EffectiveElectronicConductivity;
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
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            
            state.(elyte).T    = state.T(model.(elyte).G.mappings.cellmap);
            state.(ne).(eac).T = state.T(model.(ne).(eac).G.mappings.cellmap);
            state.(ne).(cc).T  = state.T(model.(ne).(cc).G.mappings.cellmap);
            state.(pe).(eac).T = state.T(model.(pe).(eac).G.mappings.cellmap);
            state.(pe).(cc).T  = state.T(model.(pe).(cc).G.mappings.cellmap);
            
            % Update temperature in the active materials of the electrodes.
            state.(ne).(eac) = model.(ne).(eac).updateTemperature(state.(ne).(eac));
            state.(pe).(eac) = model.(pe).(eac).updateTemperature(state.(pe).(eac));
            
        end
        
        function model = setElectrolyteVolumeFraction(model)
        % Abbreviations used in this function 
        % elyte : Electrolyte
        % ne    : NegativeElectrode
        % pe    : PositiveElectrode
            
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            eac   = 'ElectrodeActiveComponent';
            sep   = 'Separator';

            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(model.(elyte).G.mappings.cellmap) = (1 : model.(elyte).G.cells.num)';

            model.(elyte).volumeFraction = NaN(model.(elyte).G.cells.num, 1);
            model.(elyte).volumeFraction(elyte_cells(model.(ne).(eac).G.mappings.cellmap))  = model.(ne).(eac).porosity;
            model.(elyte).volumeFraction(elyte_cells(model.(pe).(eac).G.mappings.cellmap))  = model.(pe).(eac).porosity;
            model.(elyte).volumeFraction(elyte_cells(model.(elyte).(sep).G.mappings.cellmap)) = model.(elyte).(sep).porosity;


        end
        
        function model = setElectrolyteEffectiveThermalConductivity(model)
            
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            eac   = 'ElectrodeActiveComponent';

            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(model.(elyte).G.mappings.cellmap) = (1 : model.(elyte).G.cells.num)';
            elyte_cells_sep = elyte_cells(model.(sep).G.mappings.cellmap);
            model.(elyte).EffectiveThermalConductivity(elyte_cells_sep)  = model.(ne).(eac).porosity.*model.(elyte).thermalConductivity;
            model.(elyte).EffectiveThermalConductivity(elyte_cells_sep)  = model.(pe).(eac).porosity.*model.(elyte).thermalConductivity;
            
        end
        
        function initstate = setupInitialState(model)
        % Setup initial state
        %
        % Abbreviations used in this function 
        % elyte : Electrolyte
        % ne    : NegativeElectrode
        % pe    : PositiveElectrode
        % eac   : ElectrodeActiveComponent
        % cc    : CurrentCollector
            
            nc = model.G.cells.num;

            SOC = model.SOC;
            T   = model.initT;
            
            initstate.T   = T*ones(nc, 1);
            initstate.SOC = SOC*ones(nc, 1);
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            
            %% synchronize temperatures
            initstate = model.updateTemperature(initstate);
            
            %% setup initial NegativeElectrode state
            
            % shortcut
            % negAm : ActiveMaterial of the negative electrode
            
            negAm = bat.(ne).(eac).(am); 
            
            m = (1 ./ (negAm.theta100 - negAm.theta0));
            b = -m .* negAm.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* negAm.Li.cmax;
            c = c*ones(negAm.G.cells.num, 1);

            initstate.(ne).(eac).c = c;
            initstate.(ne).(eac) = model.(ne).(eac).updateChargeCarrier(initstate.(ne).(eac));
            initstate.(ne).(eac).(am) = negAm.updateOCP(initstate.(ne).(eac).(am));

            OCP = initstate.(ne).(eac).(am).OCP;
            ref = OCP(1);
            
            initstate.(ne).(eac).phi = OCP - ref;

            %% setup initial PositiveElectrode state

            % shortcut
            % posAm : ActiveMaterial of the positive electrode
            
            posAm = bat.(pe).(eac).(am);
            
            m = (1 ./ (posAm.theta100 - posAm.theta0));
            b = -m .* posAm.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* posAm.Li.cmax;
            c = c*ones(posAm.G.cells.num, 1);

            initstate.(pe).(eac).c = c;
            initstate.(pe).(eac) = model.(pe).(eac).updateChargeCarrier(initstate.(pe).(eac));
            initstate.(pe).(eac).(am) = posAm.updateOCP(initstate.(pe).(eac).(am));
            
            OCP = initstate.(pe).(eac).(am).OCP;
            initstate.(pe).(eac).phi = OCP - ref;

            %% setup initial Electrolyte state

            initstate.(elyte).phi = zeros(bat.(elyte).G.cells.num, 1)-ref;
            cs = cell(2,1);
            initstate.(elyte).cs = cs;
            initstate.(elyte).cs{1} = 1000*ones(bat.(elyte).G.cells.num, 1);

            %% setup initial Current collectors state

            OCP = initstate.(ne).(eac).(am).OCP;
            OCP = OCP(1) .* ones(bat.(ne).(cc).G.cells.num, 1);
            initstate.(ne).(cc).phi = OCP-ref;

            OCP = initstate.(pe).(eac).(am).OCP;
            OCP = OCP(1) .* ones(bat.(pe).(cc).G.cells.num, 1);
            initstate.(pe).(cc).phi = OCP - ref;
            
            initstate.(pe).(cc).E = OCP(1) - ref;
            
        end
        
        function state = updateElectrolyteCoupling(model, state)
        % Setup the electrolyte coupling by adding ion sources from the electrodes
        % shortcuts:
        % c_source : Source term for charge carrier.
                        
            battery = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            
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
            coupterm = getCoupTerm(couplingterms, 'NegativeElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = ne_R.*vols(elytecells); % we divide with F later
            
            pe_R = state.(pe).(eac).(am).R;
            coupterm = getCoupTerm(couplingterms, 'PositiveElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = pe_R.*vols(elytecells);
            
            elyte_e_source = elyte_c_source.*battery.(elyte).sp.z{1}; % we divide with F later
            
            state.Electrolyte.(ccSourceName) = elyte_c_source/F; 
            state.Electrolyte.eSource = elyte_e_source;
            
        end
        
        function state = updateAccumTerms(model, state, state0, dt)
                    
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            
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

        function state = updateElyteHeatBcSource(model, state)
        % trivial update to zero.
            elyte = 'Electrolyte';
            nc = model.(elyte).G.cells.num;

            state.(elyte).jHeatBcSource = zeros(nc, 1);
            
        end
        
        
        function state = updateEnergyAccumTerms(model, state, state0, dt)
                    
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            
            submodels = {model.(elyte)   , ...
                         model.(ne).(eac), ...
                         model.(ne).(cc) , ...
                         model.(pe).(eac), ...
                         model.(pe).(cc)};
            
            fds = {{elyte, 'accumHeat'}  , ...
                   {ne, eac, 'accumHeat'}, ...
                   {ne, cc, 'accumHeat'} , ...
                   {pe, eac, 'accumHeat'}, ...
                   {pe, cc, 'accumHeat'}};
            
            globT = state.T;
            globT0 = state0.T;
            
            for ind = 1 : numel(fds)
                
                submodel = submodels{ind};
                hC       = submodel.heatCapacity;
                if isprop(submodel, 'volumeFraction')
                    volfrac  = submodel.volumeFraction;
                else
                    volfrac = 1;
                end
                vols     = submodel.G.cells.volumes;
                cellmap  = submodel.G.mappings.cellmap;
                
                T  = globT(cellmap);
                T0 = globT0(cellmap);
                
                accumHeat = hC.*volfrac.*vols.*(T - T0)/dt;
                
                state = model.setProp(state, fds{ind}, accumHeat);
            end            
            
        end

        
        function state = updateElectrodeCoupling(model, state)
        % Setup electrode coupling by updating the potential and concentration of the electrolyte in the active component of the
        % electrode. There, those quantities are considered as input and used to compute the reaction rate.
        %
        %
        % WARNING : at the moment, we do not pass the concentrations
        %
        % shortcuts:
        % elyte : Electrolyte
        % neac  : NegativeElectrode.ElectrodeActiveComponent 
        % peac  : PositiveElectrode.ElectrodeActiveComponent
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            
            phi_elyte = state.(elyte).phi;
            
            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(bat.(elyte).G.mappings.cellmap) = (1 : bat.(elyte).G.cells.num)';

            phi_elyte_neac = phi_elyte(elyte_cells(bat.(ne).(eac).G.mappings.cellmap));
            phi_elyte_peac = phi_elyte(elyte_cells(bat.(pe).(eac).G.mappings.cellmap));

            state.(ne).(eac).(am).phiElectrolyte = phi_elyte_neac;
            state.(pe).(eac).(am).phiElectrolyte = phi_elyte_peac;
            
        end

        function state = setupExternalCouplingNegativeElectrode(model, state);
            
            ne = 'NegativeElectrode';
            cc = 'CurrentCollector';
           
            phi = state.(ne).(cc).phi;

            jExternal = model.(ne).(cc).setupExternalCoupling(phi, 0);
            
            state.(ne).(cc).jExternal = jExternal;
            
        end
        
        function state = setupExternalCouplingPositiveElectrode(model, state);
            
            pe = 'PositiveElectrode';
            cc = 'CurrentCollector';
           
            phi = state.(pe).(cc).phi;
            E = state.(pe).(cc).E;
            
            jExternal = model.(pe).(cc).setupExternalCoupling(phi, E);
            
            state.(pe).(cc).jExternal = jExternal;
            
        end

        
        function state =updateEnergyConservation(model, state);
            
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            
            % We recover energy equation for electrolyte
            
            energyCons = state.(elyte).energyCons;
            
            %  we have to add the contribution from the electrode's active component
           
            eldes = {ne, pe};
            
            for ind = 1 : numel(eldes)
                elde = eldes{ind};
                cellmap = model.mappings.(elde);
                elde_energyCons = state.(elde).(eac).energyCons;
                energyCons(cellmap) = energyCons(cellmap) + elde_energyCons;
            end
            
            % We retrieve the energy equation from the current collectors.
            
            necc_energyCons = state.(ne).(cc).energyCons;
            pecc_energyCons = state.(pe).(cc).energyCons;
            
            % We concatenate the contributions
            
            state.energyCons = vertcat(necc_energyCons, ...
                                       energyCons, ...
                                       pecc_energyCons);
        
        end
        
        function state = initStateAD(model,state)
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            
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
             state.T              , ...
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
                    state.T              , ...
                    state.(pe).(cc).E, ...
                    opts); 
            % PRIMARY variables
        end
        
        function p = getPrimaryVariables(model)
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            
            p = {{elyte, 'cs', 1} , ...
                 {elyte, 'phi'}   , ...   
                 {ne, eac, 'c'}   , ...    
                 {ne, eac, 'phi'} , ...   
                 {pe, eac, 'c'}   , ...    
                 {pe, eac, 'phi'} , ...   
                 {ne, cc, 'phi'}  , ...    
                 {pe, cc, 'phi'}  , ...
                 {'T'}            , ...
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
              {'PositiveElectrode','ElectrodeActiveComponent'}, ...
              {'NegativeElectrode','ElectrodeActiveComponent'}, ...
              {'NegativeElectrode','CurrentCollector'}, ...
              {'PositiveElectrode','CurrentCollector'}};
            model.Electrolyte.AutoDiffBackend=model.AutoDiffBackend;
            model.Electrolyte=model.Electrolyte.validateModel(varargin{:});
            model.PositiveElectrode.ElectrodeActiveComponent.AutoDiffBackend= model.AutoDiffBackend;
            model.PositiveElectrode.ElectrodeActiveComponent = model.PositiveElectrode.ElectrodeActiveComponent.validateModel(varargin{:});
            model.NegativeElectrode.ElectrodeActiveComponent.AutoDiffBackend= model.AutoDiffBackend;
            model.NegativeElectrode.ElectrodeActiveComponent = model.NegativeElectrode.ElectrodeActiveComponent.validateModel(varargin{:});
            model.NegativeElectrode.CurrentCollector.AutoDiffBackend= model.AutoDiffBackend;
            model.NegativeElectrode.CurrentCollector=model.NegativeElectrode.CurrentCollector.validateModel(varargin{:});
            model.PositiveElectrode.CurrentCollector.AutoDiffBackend=model.AutoDiffBackend;
            model.PositiveElectrode.CurrentCollector= model.PositiveElectrode.CurrentCollector.validateModel(varargin{:});
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
